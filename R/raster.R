## Raster matrix track (additive): a cells x bins matrix drawn as one raster image per window, e.g. a copy-number
## heat map of thousands of cells under the genome tracks.
##
## Data: a GRanges of bins whose metadata hold the matrix, either as one matrix column (formatting raster.field,
## rows = bins, columns = cells) or as plain numeric columns (one per cell; raster.field = NA uses them all).
## Formatting: raster = TRUE; raster.col = comma-separated palette (default blue,white,red); raster.mid = value drawn
## with the middle color (e.g. 2 copies); cmap.min / cmap.max saturate the scale (default 1% / 99% quantiles);
## raster.px = image columns across the plot (default 1500); raster.labels = draw row names (TRUE when <= 40 rows);
## raster.legend = color bar at the right (default TRUE). Row 1 is drawn at the top.

.gt_raster_matrix = function(gr, field = NA)
{
  mc = S4Vectors::mcols(gr)
  if (!is.na(field) && field %in% names(mc))
  {
    m = as.matrix(mc[[field]])
    if (is.null(colnames(m))) colnames(m) = paste0('row', seq_len(ncol(m)))
  } else {
    num = names(mc)[vapply(seq_len(ncol(mc)), function(k) is.numeric(mc[[k]]) && is.null(dim(mc[[k]])), TRUE)]
    if (!length(num)) stop('raster track: no numeric columns (or raster.field matrix) in the data')
    m = as.matrix(as.data.frame(mc[, num, drop = FALSE]))
  }
  m ## bins x rows
}

.gt_raster_colors = function(v, cmin, cmax, cols, mid = NA, na.col = 'grey92', n = 101)
{
  pal = grDevices::colorRampPalette(cols)(n)
  v = pmin(pmax(v, cmin), cmax)
  t = if (!is.na(mid) && mid > cmin && mid < cmax)
        ifelse(v <= mid, 0.5 * (v - cmin) / (mid - cmin), 0.5 + 0.5 * (v - mid) / (cmax - mid))
      else (v - cmin) / max(cmax - cmin, .Machine$double.eps)
  out = pal[1 + round(t * (n - 1))]
  out[is.na(out)] = na.col
  out
}

draw.raster = function(grl, y, windows, win.gap = NULL, raster.field = NA, raster.col = NA, raster.mid = NA,
                       cmap.min = NA, cmap.max = NA, raster.px = NA, raster.labels = NA, raster.legend = NA,
                       xlim = c(0, 20000), ...)
{
  ylim.subplot = if (is.list(y)) c(y$start[1], y$end[1]) else c(y[[1]], y[[2]])
  gr = if (inherits(grl, 'GRangesList')) grl.unlist(grl) else grl
  m = .gt_raster_matrix(gr, raster.field)
  cols = if (is.na(raster.col)) c('#2166AC', 'white', '#B2182B') else trimws(strsplit(raster.col, ',')[[1]])
  if (is.na(raster.px)) raster.px = 1500
  if (is.na(raster.legend)) raster.legend = TRUE
  nr = ncol(m)
  if (is.na(raster.labels)) raster.labels = nr <= 40
  if (is.null(win.gap)) win.gap = mean(width(windows)) * 0.2
  g0 = GenomicRanges::granges(gr); GenomicRanges::strand(g0) = '*'
  mapped = gr.flatmap(g0, windows, win.gap)
  segs = mapped$grl.segs; ws = mapped$window.segs
  winlim = range(c(ws$start, ws$end + 1))
  ws$start = affine.map(ws$start, winlim, ylim = xlim); ws$end = affine.map(ws$end + 1, winlim, ylim = xlim)
  if (nrow(segs))
  {
    segs$pos1 = affine.map(segs$pos1, winlim, ylim = xlim); segs$pos2 = affine.map(segs$pos2 + 1, winlim, ylim = xlim)
  }
  vals = m[unique(segs$query.id), , drop = FALSE]
  if (is.na(cmap.min)) cmap.min = stats::quantile(vals, 0.01, na.rm = TRUE)
  if (is.na(cmap.max)) cmap.max = stats::quantile(vals, 0.99, na.rm = TRUE)
  if (!is.finite(cmap.min) || !is.finite(cmap.max)) { cmap.min = 0; cmap.max = 1 }
  if (cmap.max <= cmap.min) cmap.max = cmap.min + 1
  tot = sum(ws$end - ws$start)
  for (k in seq_len(nrow(ws)))
  {
    s = segs[segs$window == k, , drop = FALSE]
    if (!nrow(s)) next
    s = s[order(s$pos1), , drop = FALSE]
    npx = max(2, round(raster.px * (ws$end[k] - ws$start[k]) / tot))
    xc = ws$start[k] + (seq_len(npx) - 0.5) / npx * (ws$end[k] - ws$start[k])
    b = findInterval(xc, s$pos1)
    b[b == 0 | xc > s$pos2[pmax(b, 1)]] = NA ## pixels in gaps between bins
    v = m[s$query.id[b], , drop = FALSE] ## npx x rows
    img = matrix(.gt_raster_colors(as.vector(v), cmap.min, cmap.max, cols, raster.mid), nrow = npx)
    img[is.na(b), ] = 'white'
    graphics::rasterImage(grDevices::as.raster(t(img)), ws$start[k], ylim.subplot[1], ws$end[k], ylim.subplot[2], interpolate = FALSE)
  }
  graphics::rect(min(ws$start), ylim.subplot[1], max(ws$end), ylim.subplot[2], border = 'grey40', lwd = 0.4)
  cex = if (!is.null(.gts())) .pt2cex(.gts()$label.size) else 0.6
  if (isTRUE(as.logical(raster.labels)))
  {
    yc = ylim.subplot[2] - (seq_len(nr) - 0.5) / nr * diff(ylim.subplot)
    graphics::text(min(ws$start) - diff(graphics::par('usr')[1:2]) * 0.005, yc, colnames(m), adj = c(1, 0.5), cex = cex)
  }
  if (isTRUE(as.logical(raster.legend)))
  {
    lv = seq(cmap.min, cmap.max, length.out = 50)
    lc = .gt_raster_colors(lv, cmap.min, cmap.max, cols, raster.mid)
    xr = max(ws$end) + diff(graphics::par('usr')[1:2]) * c(0.01, 0.02)
    yb = ylim.subplot[1] + diff(ylim.subplot) * c(0.1, 0.9)
    graphics::rasterImage(grDevices::as.raster(matrix(rev(lc), ncol = 1)), xr[1], yb[1], xr[2], yb[2], interpolate = FALSE)
    lab = c(cmap.min, if (!is.na(raster.mid) && raster.mid > cmap.min && raster.mid < cmap.max) raster.mid, cmap.max)
    graphics::text(xr[2] + diff(graphics::par('usr')[1:2]) * 0.004, yb[1] + (lab - cmap.min) / (cmap.max - cmap.min) * diff(yb),
                   signif(lab, 3), adj = c(0, 0.5), cex = cex)
  }
  ws
}
