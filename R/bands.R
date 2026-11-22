## Row bands (additive): a colored strip and a label at the left of each track row (e.g. clade colors), and brackets
## that group several rows (e.g. the cells of one clade under rep = "cells"). Driven by formatting fields:
##   band.col   strip color per track (NA = no strip)
##   band.label text left of the strip (NA = none)
##   band.group group id; consecutive tracks with the same id get one bracket, labeled with the id
##   band.group.col bracket / group label color (first track of the group)
## Drawn after all tracks from the layout plot.gTrack returns; nothing is drawn when no band.* field is set.
.gt_draw_bands = function(.Object, lay, band.params = list())
{
  f = formatting(.Object)
  fld = function(k) if (is.null(f[[k]])) rep(NA, length(.Object)) else f[[k]]
  bcol = as.character(fld('band.col')); blab = as.character(fld('band.label')); bgrp = as.character(fld('band.group'))
  gcol = as.character(fld('band.group.col'))
  if (all(is.na(bcol)) && all(is.na(blab)) && all(is.na(bgrp))) return(invisible(NULL))
  tr = lay$tracks
  if (!nrow(tr)) return(invisible(NULL))
  P = utils::modifyList(list(offset.mm = 9, strip.mm = 1.2, bracket.mm = 1.2, label.mm = NA, cex = NA), band.params)
  st = .gts()
  cex = if (!is.na(P$cex)) P$cex else if (!is.null(st)) .pt2cex(st$label.size) else 0.7
  x0 = min(lay$windows$x0)
  xr = x0 - .mm2x(P$offset.mm) ## right edge of the strip
  xs = xr - .mm2x(P$strip.mm)
  has.strip = any(!is.na(bcol))
  for (k in seq_len(nrow(tr)))
  {
    j = tr$track[k]
    if (!is.na(bcol[j])) graphics::rect(xs, tr$y0[k], xr, tr$y1[k], col = bcol[j], border = NA)
    if (!is.na(blab[j]))
      graphics::text((if (has.strip) xs else xr) - .mm2x(0.6), (tr$y0[k] + tr$y1[k]) / 2, blab[j], adj = c(1, 0.5), cex = cex,
                     col = if (!is.na(bcol[j]) && is.na(bgrp[j])) bcol[j] else 'black')
  }
  if (any(!is.na(bgrp)))
  {
    lw = if (!is.na(P$label.mm)) .mm2x(P$label.mm) else if (any(!is.na(blab))) max(graphics::strwidth(blab[!is.na(blab)], cex = cex)) + .mm2x(0.6) else 0
    xb = (if (has.strip) xs else xr) - lw - .mm2x(P$bracket.mm)
    r = rle(ifelse(is.na(bgrp[tr$track]), '\001', bgrp[tr$track])); ends = cumsum(r$lengths); starts = ends - r$lengths + 1
    for (g in seq_along(r$values))
    {
      if (r$values[g] == '\001') next
      ix = starts[g]:ends[g]
      y0 = min(tr$y0[ix]); y1 = max(tr$y1[ix]); cl = gcol[tr$track[ix[1]]]; if (is.na(cl)) cl = 'black'
      tk = .mm2x(0.8)
      graphics::segments(c(xb, xb, xb), c(y0, y0, y1), c(xb, xb + tk, xb + tk), c(y1, y0, y1), col = cl, lwd = 0.8, lend = 1)
      graphics::text(xb - .mm2x(0.8), (y0 + y1) / 2, r$values[g], srt = 90, adj = c(0.5, 0), cex = cex, col = cl)
    }
  }
  invisible(NULL)
}
