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
#' @param win.widths relative display widths of the windows (NULL = proportional to genomic span, the classic layout);
#'   when given (one value per window), every window is drawn at its share of the plot width, so a small amplicon window is
#'   not a sliver; data inside a window keep their linear genomic scale; edges and paths between windows are drawn as usual
#' @param xaxis.chr.merge write the chromosome label once under each run of consecutive windows on the same chromosome
#' @param arc.style "gtrack" (classic connectors) or "bell": every ALT edge (edges with $not.flat) is one bell that rises
#'   vertically from both segment ends to an apex height proportional to log10 of the junction span (edges$span, bp; else the
#'   x distance), clipped to the top of the track; REF edges are straight thin connectors between the segment ends; edges to
#'   segments outside the windows are a short arrow at the segment end labelled edges$dangle.label (e.g. "chr2:33 Mb")
#' @param arc.band apex height of the tallest bell as a fraction of the track's y range
#' @param arc.span spans (bp) that map to the lowest and the highest bell
#' @param arc.range bell heights (fractions of arc.band) for the spans in arc.span
#' @param label.halo white halo behind edge labels; with arc.style "bell" edge labels sit above the bell apex and are moved
#'   apart (simple repel) when they would overlap
#' @param dangle.mm length of the dangling-edge arrow (x, y) in mm
#' @return list
#' @export
pub_style = function(family = "Helvetica", fontsize = 7, title.size = 8, title.font = 1, label.size = 6,
                     lwd.axis = 0.5, tick.mm = 1, ylab.gap.mm = 0.6,
                     xaxis.unit = 1e6, xaxis.suffix = "Mb", xaxis.nticks = 5, xaxis.chr.prefix = "chr",
                     yaxis.nticks = 2, sep.col = "grey45", sep.lty = 2, sep.lwd = 0.4,
                     mar.mm = c(6.5, 7, 1, 2), ygrid = FALSE, win.widths = NULL, xaxis.chr.merge = FALSE,
                     arc.style = "gtrack", arc.band = 0.35, arc.span = c(1e3, 1e8), arc.range = c(0.3, 1), label.halo = TRUE,
                     dangle.mm = c(1.6, 1.4))
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
    chrs = as.character(GenomicRanges::seqnames(windows))
    xc = (x0 + x1) / 2
    if (isTRUE(st$xaxis.chr.merge)) ## one label per run of consecutive windows on the same chromosome
    {
      run = cumsum(c(TRUE, chrs[-1] != chrs[-nwin]))
      if (k > 1 && run[k] == run[k - 1]) next
      rk = which(run == run[k])
      xc = (window.segs$start[min(rk)] + window.segs$end[max(rk)]) / 2
    }
    if (!is.na(st$xaxis.chr.prefix) && !grepl("^chr", chr))
      chr = paste0(st$xaxis.chr.prefix, chr)
    if (!is.na(st$xaxis.suffix))
      chr = paste0(chr, " (", st$xaxis.suffix, ")")
    graphics::text(xc, xaxis.pos - tl - .mm2y(0.5) - 1.15 * graphics::strheight("0", cex = cex) - .mm2y(0.6), chr,
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
  lab = names(line.loc)
  if (length(line.loc) > 2) ## thin labels that would overlap (keep the first and every k-th)
  {
    h = 1.15 * graphics::strheight("0", cex = cex)
    k = ceiling(h / min(diff(sort(line.loc))))
    if (k > 1)
      lab[(seq_along(lab) - 1) %% k != 0] = ""
  }
  graphics::text(x - tl - .mm2x(st$ylab.gap.mm), line.loc, lab, adj = c(1, 0.5), cex = cex)
}

## internal: text with a white halo (8 offset copies underneath)
.halo_text = function(x, y, labels, col, cex, adj, halo = "white", r.mm = 0.22)
{
  if (!length(x)) return(invisible())
  for (th in seq(0, 2 * pi, length.out = 9)[-9])
    graphics::text(x + .mm2x(r.mm) * cos(th), y + .mm2y(r.mm) * sin(th), labels, col = halo, cex = cex, adj = adj)
  graphics::text(x, y, labels, col = col, cex = cex, adj = adj)
}

## internal: place labels (box centre x, bottom y) so that they overlap neither each other nor the lines in `obst` (matrix of x, y points
## sampled along the drawn arcs): try the spot at the anchor, then beside it, then higher; boxes are kept inside [xr] and below ymax
.repel = function(x, y, w, h, ymax = Inf, xr = c(-Inf, Inf), obst = NULL, vert = rep(FALSE, length(x)))
{
  n = length(x); px = x; py = y; placed = NULL
  if (!is.matrix(xr)) xr = matrix(xr, nrow = n, ncol = 2, byrow = TRUE)  ## per-label x limits
  clampx = function(cx, k) min(max(cx, xr[k, 1] + w[k] / 2), xr[k, 2] - w[k] / 2)
  ## cost of a box: overlapping labels count 100, every arc point inside the box 1 (the box is padded by 15 % of its height)
  cost = function(cx, cy, k)
    (if (is.null(placed)) 0 else 100 * sum(abs(placed[, 1] - cx) < (placed[, 3] + w[k]) / 2 & abs(placed[, 2] - cy) < (placed[, 4] + h[k]) / 2)) +
    (if (is.null(obst)) 0 else sum(abs(obst[, 1] - cx) < w[k] / 2 & obst[, 2] > cy - 0.15 * h[k] & obst[, 2] < cy + 1.15 * h[k]))
  cand = rbind(c(0, 0), c(0, 0.6), c(0.6, 0), c(-0.6, 0), c(0, 1.2), c(0.6, 0.8), c(-0.6, 0.8), c(1.1, 0), c(-1.1, 0), c(0, 1.8), c(1.1, 1), c(-1.1, 1), c(0, 2.6))
  cand.v = rbind(cbind(0, 0:5), cbind(0.5, 0:3), cbind(-0.5, 0:3))  ## dangling-edge labels: stack upwards first
  for (k in order(vert, -y))
  {
    best = NULL; bc = Inf; cc = if (vert[k]) cand.v else cand
    for (m in seq_len(nrow(cc)))
    {
      cx = clampx(x[k] + cc[m, 1] * w[k], k); cy = min(y[k] + cc[m, 2] * h[k], ymax - h[k])
      co = cost(cx, cy, k) + 0.01 * m  ## ties: the earlier (closer) candidate
      if (co < bc) { bc = co; best = c(cx, cy) }
      if (co < 1) break
    }
    px[k] = best[1]; py[k] = best[2]; placed = rbind(placed, c(best[1], best[2], w[k], h[k]))
  }
  list(x = px, y = py)
}

## internal: points along the 6-point bezier that connectors() draws for a bell (type "U", h = 0) from (x0, y0) to (x1, y1) with top control y = yt
.bell_pts = function(x0, y0, x1, y1, yt, n = 40)
{
  t = seq(0, 1, length.out = n); B = sapply(0:5, function(i) choose(5, i) * t^i * (1 - t)^(5 - i))
  do.call(rbind, lapply(seq_along(x0), function(k) cbind(B %*% c(x0[k], x0[k], x0[k], x1[k], x1[k], x1[k]), B %*% c(y0[k], y0[k], yt[k], yt[k], y1[k], y1[k]))))
}

## internal: edges in the "bell" style (pub_style(arc.style = "bell")), called by draw.grl instead of the classic connectors
.pub_draw_edges = function(st, edges, yr, ws = NULL)
{
  alt = if (is.null(edges$not.flat)) rep(FALSE, nrow(edges)) else edges$not.flat %in% TRUE
  dg = edges$dangle %in% TRUE
  cexl = .pt2cex(st$label.size)
  lab.h = graphics::strheight("J", cex = cexl)
  lb = if (is.null(edges$label)) rep(NA_character_, nrow(edges)) else as.character(edges$label)
  ## REF: straight thin connectors between the segment ends (vertical between adjacent segments); dangling REF edges are not drawn
  r = which(!alt & !dg)
  if (length(r))
    graphics::segments(edges$x.pos.from[r], edges$y.pos.from[r], edges$x.pos.to[r], edges$y.pos.to[r],
                       col = edges$col[r], lwd = edges$lwd[r], lty = edges$lty[r], lend = 1)
  ## ALT: one bell per junction, from the two segment ends straight up to an apex at height ~ log10(span), kept inside the track
  a = which(alt & !dg)
  apex = rep(NA_real_, nrow(edges))
  if (length(a))
  {
    band = st$arc.band * diff(yr)
    sp = if (!is.null(edges$span)) suppressWarnings(as.numeric(edges$span[a])) else rep(NA_real_, length(a))
    fx = abs(edges$x.pos.to[a] - edges$x.pos.from[a]) / diff(graphics::par("usr")[1:2])
    f = ifelse(is.na(sp), affine.map(fx, xlim = c(0, 1), ylim = st$arc.range, cap = TRUE),
               affine.map(log10(pmax(sp, 1)), xlim = log10(st$arc.span), ylim = st$arc.range, cap = TRUE))
    base = pmax(edges$y.pos.from[a], edges$y.pos.to[a])
    room = yr[2] - base - 1.5 * lab.h
    h = pmax(pmin(band * f, room), 0.12 * band)
    connectors(edges$x.pos.from[a], edges$y.pos.from[a], edges$dir.from[a], edges$x.pos.to[a], edges$y.pos.to[a], edges$dir.to[a],
               v = h / 0.625, h = 0, type = "U", f.arrow = FALSE, b.arrow = FALSE, cex.arrow = 0, col.arrow = edges$col[a],
               lwd = edges$lwd[a], lty = edges$lty[a], col = edges$col[a])
    apex[a] = base + h
    obst = .bell_pts(edges$x.pos.from[a], edges$y.pos.from[a], edges$x.pos.to[a], edges$y.pos.to[a], base + h / 0.625)
  }
  else obst = NULL
  ## dangling ALT edges (mate segment outside the windows): short arrow from the segment end, pointing out and down; the mate locus
  ## (edges$dangle.label) is written in a row under the track (up to 3 rows, inside the arrow's window), with a thin leader to the arrow tip
  d = which(alt & dg)
  if (length(d))
  {
    xa = ifelse(edges$from.ix[d], edges$x.pos.from[d], edges$x.pos.to[d]); ya = ifelse(edges$from.ix[d], edges$y.pos.from[d], edges$y.pos.to[d])
    side = ifelse(edges$from.ix[d], edges$dir.from[d], edges$dir.to[d]); dir = ifelse(side == '+', 1, -1)  ## right end of a segment -> points right
    dlab = if (is.null(edges$dangle.label)) rep(NA_character_, length(d)) else as.character(edges$dangle.label[d])
    u = !duplicated(data.frame(round(xa, 3), round(ya, 5), dir, dlab))  ## gGraph tracks hold both strand copies of every edge
    d = d[u]; xa = xa[u]; ya = ya[u]; dir = dir[u]; dlab = dlab[u]
    x1 = xa + dir * .mm2x(st$dangle.mm[1]); y1 = ya - .mm2y(st$dangle.mm[2])
    graphics::arrows(xa, ya, x1, y1, length = 0.03, angle = 25, col = edges$col[d], lwd = 0.7 * edges$lwd[d])
    dl = data.frame(x = x1, y = y1, lab = dlab, col = as.character(edges$col[d]), stringsAsFactors = FALSE)
    dl = dl[!is.na(dl$lab) & nzchar(dl$lab), , drop = FALSE]
    dl = dl[order(dl$x), , drop = FALSE]
    if (nrow(dl) > 1)  ## one label for arrows to the same mate within 4 mm (eg both ends of an inserted fragment)
    {
      keep = rep(TRUE, nrow(dl))
      for (i in seq_len(nrow(dl))[-1]) if (any(keep[seq_len(i - 1)] & dl$lab[seq_len(i - 1)] == dl$lab[i] & abs(dl$x[seq_len(i - 1)] - dl$x[i]) < .mm2x(4))) keep[i] = FALSE
      dl = dl[keep, , drop = FALSE]
    }
    if (nrow(dl))
    {
      w = graphics::strwidth(dl$lab, cex = cexl) * 1.06; gap = .mm2x(0.8); rh = lab.h * 1.2
      u = graphics::par("usr")[1:2]; X = u + c(1, -1) * diff(u) * 0.04 / 1.08
      wk = if (!is.null(ws) && nrow(ws)) vapply(dl$x, function(x) { k = which(ws$start <= x & ws$end >= x)[1]; if (is.na(k)) which.min(abs((ws$start + ws$end) / 2 - x)) else k }, 1L) else rep(1L, nrow(dl))
      re = list()  ## right edge of the last label per window and row
      for (i in seq_len(nrow(dl)))
      {
        xr = if (!is.null(ws) && nrow(ws)) c(ws$start[wk[i]], ws$end[wk[i]]) else X
        for (r in 1:3)
        {
          key = paste(wk[i], r); last = if (is.null(re[[key]])) -Inf else re[[key]]
          cx = min(max(dl$x[i], last + gap + w[i] / 2, xr[1] + w[i] / 2), xr[2] - w[i] / 2)
          if (cx - w[i] / 2 >= last + gap - 1e-9 || r == 3) break
        }
        re[[key]] = cx + w[i] / 2
        ytop = yr[1] - .mm2y(0.4) - (r - 1) * rh
        if (abs(cx - dl$x[i]) > .mm2x(0.6) || r > 1)
          graphics::segments(dl$x[i], dl$y[i], cx, ytop, col = dl$col[i], lwd = 0.4)
        if (isTRUE(st$label.halo)) .halo_text(cx, ytop, dl$lab[i], col = dl$col[i], cex = cexl, adj = c(0.5, 1))
        else graphics::text(cx, ytop, dl$lab[i], col = dl$col[i], cex = cexl, adj = c(0.5, 1))
      }
    }
  }
  ## edge labels above the bell apex; white halo; moved apart when they overlap each other or an arc
  l = which(!is.na(apex) & !is.na(lb) & nzchar(lb))
  lab.col = if (!is.null(edges$label.col)) edges$label.col else edges$col
  L = data.frame(x = (edges$x.pos.from[l] + edges$x.pos.to[l]) / 2, y = apex[l] + .mm2y(0.35), lab = lb[l],
                 col = as.character(lab.col[l]), dir = rep(0, length(l)), stringsAsFactors = FALSE)
  L = L[!duplicated(L[, c("x", "y", "lab")]), , drop = FALSE]  ## strand copies of one edge
  if (nrow(L))
  {
    w = graphics::strwidth(L$lab, cex = cexl) * 1.08; hh = rep(lab.h * 1.15, nrow(L))
    cx = L$x + ifelse(L$dir < 0, -w / 2, ifelse(L$dir > 0, w / 2, 0))  ## centre of each label box
    u = graphics::par("usr")[1:2]; xr = matrix(u + c(1, -1) * diff(u) * 0.04 / 1.08, nrow = nrow(L), ncol = 2, byrow = TRUE)  ## the windows (plot() pads x by 4 %)
    if (!is.null(ws) && nrow(ws)) for (i in which(L$dir != 0)) {  ## a dangling-edge label stays inside the window of its arrow
      k = which(ws$start <= L$x[i] - L$dir[i] * .mm2x(0.5) & ws$end >= L$x[i] - L$dir[i] * .mm2x(0.5))[1]
      if (!is.na(k)) xr[i, ] = c(ws$start[k], ws$end[k]) }
    p = .repel(cx, L$y, w, hh, ymax = yr[2] + lab.h, xr = xr, obst = obst, vert = L$dir != 0)
    if (isTRUE(st$label.halo))
      .halo_text(p$x, p$y, L$lab, col = L$col, cex = cexl, adj = c(0.5, 0))
    else
      graphics::text(p$x, p$y, L$lab, col = L$col, cex = cexl, adj = c(0.5, 0))
  }
  invisible(apex)
}

#' @name pub_window_layout
#' @title Where plot.gTrack puts each window
#' @description
#' Horizontal position of every window of plot(gt, windows, style = style) as a fraction of the plot region width (the device
#' width minus style$mar.mm[2] and style$mar.mm[4]); honours style$win.widths.  Used to place overlays or sidecar hit regions.
#' @param windows GRanges windows as given to plot.gTrack
#' @param style pub_style() list
#' @param gap gap between windows in bp (plot.gTrack default: total window width / 30)
#' @return data.frame chr, start, end, x0, x1
#' @export
pub_window_layout = function(windows, style = pub_style(), gap = NULL)
{
  old = .gt_style$style; .gt_style$style = style; on.exit(.gt_style$style <- old, add = TRUE)
  if (is.null(gap)) gap = sum(as.numeric(GenomicRanges::width(windows))) / 30
  w = windows; suppressWarnings(GenomicRanges::end(w) <- GenomicRanges::end(w) + 1)
  ws = gr.flatmap(GenomicRanges::GRanges(), w, gap)$window.segs
  lim = range(c(ws$start, ws$end)); pad = 0.04  ## plot() default xaxs = "r"
  data.frame(chr = as.character(GenomicRanges::seqnames(windows)), start = GenomicRanges::start(windows), end = GenomicRanges::end(windows),
             x0 = ((ws$start - lim[1]) / diff(lim) + pad) / (1 + 2 * pad), x1 = ((ws$end - lim[1]) / diff(lim) + pad) / (1 + 2 * pad))
}
