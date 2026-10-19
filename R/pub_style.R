## Publication style for gTrack plots.
##
## pub_style() returns a list of style settings; pass it to plot(gt, win, style = pub_style()) or set it
## globally with options(gTrack.style = pub_style()).  When a style is active, plot.gTrack
##   * uses the given font family and absolute font sizes (pt) instead of cex multiples,
##   * drops the grey window background and the y grid, draws thin y axes with outward ticks of a fixed
##     length (mm) and horizontal labels,
##   * draws one clean x axis at the bottom with "pretty" ticks in Mb and the chromosome name below,
##   * writes track titles left-aligned above each track,
##   * colours junction (edge) labels like their edge and places them at the arc apex.
## Without a style (the default) gTrack draws exactly as before.
## gtrack_grob() renders a gTrack plot as a vector grid grob (via gridGraphics) for use in grid layouts.

.gt_style = new.env()

#' @name pub_style
#' @title Publication style settings for gTrack
#' @description
#' List of style settings consumed by plot.gTrack(style = ...).
#' @param family font family (e.g. "Helvetica" for pdf(), "Nimbus Sans" or "Helvetica" for ragg PNG)
#' @param fontsize base font size in pt (axis tick labels)
#' @param title.size track title size in pt
#' @param title.font track title font (1 plain, 2 bold, 3 italic)
#' @param label.size size in pt of edge and range labels
#' @param lwd.axis axis line width in pt
#' @param tick.mm tick length in mm (ticks point outwards)
#' @param ylab.gap.mm gap between y axis line and tick labels in mm
#' @param xaxis.unit genomic unit of the x axis (1e6 = Mb)
#' @param xaxis.suffix unit appended to the chromosome label under each window, eg "chr7 (Mb)" (NA for none)
#' @param xaxis.nticks approximate number of x ticks per window
#' @param xaxis.chr.prefix prefix for chromosome labels under each window ("chr" gives "chr7")
#' @param yaxis.nticks approximate number of y ticks
#' @param sep.col,sep.lty,sep.lwd window separator line style (drawn between windows only)
#' @param mar.mm plot margins (bottom, left, top, right) in mm (NULL keeps the caller's par("mar"))
#' @param ygrid draw horizontal y grid lines (FALSE)
#' @return list
#' @export
pub_style = function(family = "Helvetica", fontsize = 7, title.size = 8, title.font = 1, label.size = 6,
                     lwd.axis = 0.5, tick.mm = 1, ylab.gap.mm = 0.6,
                     xaxis.unit = 1e6, xaxis.suffix = "Mb", xaxis.nticks = 5, xaxis.chr.prefix = "chr",
                     yaxis.nticks = 2, sep.col = "grey45", sep.lty = 2, sep.lwd = 0.4,
                     mar.mm = c(6.5, 7, 1, 2), ygrid = FALSE)
{
  structure(as.list(environment()), class = c("gTrackStyle", "list"))
}

## internal: currently active style (NULL if none)
.gts = function() .gt_style$style

## internal: font size in pt -> cex relative to current par('ps')
.pt2cex = function(pt) pt / graphics::par("ps")

## internal: multiplier for cex-specified labels (cex = 1 means style$fontsize pt in pub style)
.cexscale = function() if (is.null(.gts())) 1 else .pt2cex(.gts()$fontsize)

## internal: mm -> user units along y / x
.mm2y = function(mm) abs(diff(graphics::grconvertY(c(0, mm / 25.4), "inches", "user")))
.mm2x = function(mm) abs(diff(graphics::grconvertX(c(0, mm / 25.4), "inches", "user")))

#' @name gtrack_grob
#' @title Render a gTrack plot as a grid grob
#' @description
#' Draws plot(gt, windows, ...) with base graphics and converts it into a grid grob (gridGraphics::echoGrob),
#' so that it can be placed into grid / patchwork layouts and stays vector in pdf output.  The grob is
#' re-drawn at the size of the viewport it is placed in.
#' @param gt gTrack object
#' @param windows windows passed to plot.gTrack
#' @param style style list (default pub_style())
#' @param ... further arguments to plot.gTrack (eg links, y.heights)
#' @return grob
#' @export
gtrack_grob = function(gt, windows, style = pub_style(), ...)
{
  if (!requireNamespace("gridGraphics", quietly = TRUE))
    stop("gtrack_grob needs the gridGraphics package")
  args = list(...)
  gridGraphics::echoGrob(function() do.call(plot.gTrack, c(list(x = gt, windows = windows, style = style), args)), name = "gtrack")
}

## internal: x axis in publication style (one axis per window, pretty ticks, chromosome name underneath)
.pub_xaxis = function(st, windows, window.segs, xaxis.pos, sep.draw = TRUE)
{
  u = st$xaxis.unit
  tl = .mm2y(st$tick.mm)
  cex = .pt2cex(st$fontsize)
  nwin = length(windows)
  for (k in seq_len(nwin))
  {
    s = GenomicRanges::start(windows)[k]; e = GenomicRanges::end(windows)[k]
    x0 = window.segs$start[k]; x1 = window.segs$end[k]
    at.g = pretty(c(s, e) / u, n = max(2, round(st$xaxis.nticks * (x1 - x0) / diff(graphics::par("usr")[1:2]))))
    at.g = at.g[at.g * u >= s & at.g * u <= e]
    if (nwin > 1) ## side-by-side windows: keep tick labels away from the window edges so neighbours do not collide
      at.g = at.g[(at.g * u - s) / (e - s) > 0.08 & (e - at.g * u) / (e - s) > 0.08]
    at.x = x0 + (at.g * u - s) / (e - s) * (x1 - x0)
    graphics::segments(x0, xaxis.pos, x1, xaxis.pos, lwd = st$lwd.axis, lend = 1)
    graphics::segments(at.x, xaxis.pos, at.x, xaxis.pos - tl, lwd = st$lwd.axis, lend = 1)
    graphics::text(at.x, xaxis.pos - tl - .mm2y(0.5), format(at.g, trim = TRUE, drop0trailing = TRUE), adj = c(0.5, 1), cex = cex)
    chr = as.character(GenomicRanges::seqnames(windows))[k]
    if (!is.na(st$xaxis.chr.prefix) && !grepl("^chr", chr))
      chr = paste0(st$xaxis.chr.prefix, chr)
    if (!is.na(st$xaxis.suffix))
      chr = paste0(chr, " (", st$xaxis.suffix, ")")
    graphics::text((x0 + x1) / 2, xaxis.pos - tl - .mm2y(0.5) - 1.15 * graphics::strheight("0", cex = cex) - .mm2y(0.6), chr,
                   adj = c(0.5, 1), cex = .pt2cex(st$title.size))
  }
}

## internal: y axis in publication style
.pub_yaxis = function(st, line.loc, x)
{
  cex = .pt2cex(st$fontsize)
  tl = .mm2x(st$tick.mm)
  graphics::segments(x, min(line.loc), x, max(line.loc), lwd = st$lwd.axis, lend = 1)
  graphics::segments(x, line.loc, x - tl, line.loc, lwd = st$lwd.axis, lend = 1)
  graphics::text(x - tl - .mm2x(st$ylab.gap.mm), line.loc, names(line.loc), adj = c(1, 0.5), cex = cex)
}
