## Layout return value of plot.gTrack (additive; nothing is drawn differently).
##
## plot.gTrack() invisibly returns a "gTrackLayout": the plotted position of every track row and of every genome
## window, in user coordinates and in normalized device coordinates (NDC, 0-1 over the whole device). Callers use
## it to align panels drawn outside gTrack (a tree at the left, ggplot strips below) to the gTrack rows and x axis.

## one row per track: user y span, y data range (when the track has a y axis), and its type
.gt_layout_row = function(.Object, j, ylim.subplot, y.ticks, tname)
{
  f = formatting(.Object)
  fj = function(k, default = NA) if (!is.null(f[[k]]) && length(f[[k]]) >= j) f[[k]][j] else default
  yr = if (length(y.ticks)) range(y.ticks) else c(NA_real_, NA_real_)
  type = if (isTRUE(as.logical(fj('triangle', FALSE)))) 'triangle'
         else if (isTRUE(as.logical(fj('raster', FALSE)))) 'raster'
         else if (!is.na(fj('y.field'))) 'y' else 'ranges'
  data.frame(track = j, name = if (is.na(tname)) NA_character_ else as.character(tname), type = type,
             y0 = ylim.subplot[j, 'start'], y1 = ylim.subplot[j, 'end'],
             ydata0 = yr[1], ydata1 = yr[2],
             y.log1p = isTRUE(as.logical(fj('y.log1p', FALSE))), log = isTRUE(as.logical(fj('log', FALSE))),
             stringsAsFactors = FALSE)
}

.gt_layout = function(rows, windows, window.segs, ylim.stacks)
{
  rows = rows[!vapply(rows, is.null, TRUE)]
  tracks = if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(tracks))
  {
    tracks$ndc.y0 = graphics::grconvertY(tracks$y0, 'user', 'ndc')
    tracks$ndc.y1 = graphics::grconvertY(tracks$y1, 'user', 'ndc')
  }
  win = data.frame(seqnames = as.character(GenomicRanges::seqnames(windows)),
                   start = GenomicRanges::start(windows), end = GenomicRanges::end(windows),
                   x0 = window.segs$start, x1 = window.segs$end, stringsAsFactors = FALSE)
  win$ndc.x0 = graphics::grconvertX(win$x0, 'user', 'ndc')
  win$ndc.x1 = graphics::grconvertX(win$x1, 'user', 'ndc')
  structure(list(tracks = tracks, windows = win, usr = graphics::par('usr'), plt = graphics::par('plt'),
                 fig = graphics::par('fig'), din = graphics::par('din'),
                 xaxis.y = ylim.stacks$xaxis.pos[1], ndc.xaxis.y = graphics::grconvertY(ylim.stacks$xaxis.pos[1], 'user', 'ndc')),
            class = 'gTrackLayout')
}

#' @name gt_layout_x
#' @title Map genome positions to plotted x
#' @description
#' Maps genomic positions to the x coordinate they were drawn at by \code{plot.gTrack}, using the layout that
#' \code{plot.gTrack} returns invisibly. Positions outside every plotted window get NA.
#' @param layout the value returned by \code{plot.gTrack} (class gTrackLayout)
#' @param seqnames character vector of chromosome names (recycled), or a GRanges (its start is used)
#' @param pos numeric positions (ignored when seqnames is a GRanges)
#' @param units 'user' (plot coordinates) or 'ndc' (0-1 across the device)
#' @return numeric vector of x positions
#' @export
gt_layout_x = function(layout, seqnames, pos = NULL, units = c('user', 'ndc'))
{
  units = match.arg(units)
  if (inherits(seqnames, 'GRanges'))
  {
    pos = GenomicRanges::start(seqnames)
    seqnames = as.character(GenomicRanges::seqnames(seqnames))
  }
  n = max(length(seqnames), length(pos))
  seqnames = rep_len(sub('^chr', '', as.character(seqnames)), n); pos = rep_len(pos, n)
  w = layout$windows
  wsq = sub('^chr', '', w$seqnames)
  out = rep(NA_real_, n)
  for (k in seq_len(nrow(w)))
  {
    ix = is.na(out) & seqnames == wsq[k] & pos >= w$start[k] & pos <= w$end[k]
    if (any(ix))
    {
      a = if (units == 'user') c(w$x0[k], w$x1[k]) else c(w$ndc.x0[k], w$ndc.x1[k])
      out[ix] = a[1] + (pos[ix] - w$start[k]) / (w$end[k] - w$start[k] + 1) * (a[2] - a[1])
    }
  }
  out
}

#' @name gt_layout_y
#' @title Map a track's data values to plotted y
#' @description
#' Maps data values of track \code{track} (index or name) to the y they were drawn at, honoring \code{y.log1p}
#' and \code{log}. Tracks without a y axis map 0-1 onto their row span.
#' @param layout the value returned by \code{plot.gTrack}
#' @param track track index or name
#' @param y numeric data values
#' @param units 'user' or 'ndc'
#' @export
gt_layout_y = function(layout, track, y, units = c('user', 'ndc'))
{
  units = match.arg(units)
  tr = layout$tracks
  r = if (is.character(track)) tr[match(track, tr$name), ] else tr[match(track, tr$track), ]
  if (!nrow(r) || is.na(r$track)) stop('track not found in layout')
  yy = y
  if (r$y.log1p) yy = log10(1 + pmax(yy, 0)) else if (r$log) yy = log10(yy)
  d = if (is.na(r$ydata0)) c(0, 1) else c(r$ydata0, r$ydata1)
  a = if (units == 'user') c(r$y0, r$y1) else c(r$ndc.y0, r$ndc.y1)
  a[1] + (yy - d[1]) / (d[2] - d[1]) * (a[2] - a[1])
}

#' @export
print.gTrackLayout = function(x, ...)
{
  cat('gTrackLayout:', nrow(x$tracks), 'tracks,', nrow(x$windows), 'windows\n')
  print(x$tracks[, c('track', 'name', 'type', 'y0', 'y1', 'ndc.y0', 'ndc.y1')], row.names = FALSE)
  invisible(x)
}
