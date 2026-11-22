## Junction annotation (additive): a small bar chart at the apex of each junction arc, e.g. the fraction of cells with
## split-read support in each clade, drawn on the genome itself instead of in a separate strip.
##
## edges$bars: list column, one numeric vector (values 0-1) per edge (NULL / NA = no bars). edges$bars.col: list column
## or character (comma-separated colors per edge, recycled over the bars); default grey40. Bar width / maximum height in
## mm: edges$bars.w.mm (0.7) and edges$bars.h.mm (2.5). Each bar sits in a light frame of full height, so 0 and 1 read.
## Works with the classic connectors and with pub_style(arc.style = "bell") (edge labels move above the bars).
## A count or other text per arc stays edges$label; arc width stays edges$lwd.

.gt_bars_h = function(edges, k)
{
  if (is.null(edges$bars)) return(rep(0, length(k)))
  h = if (!is.null(edges$bars.h.mm)) as.numeric(edges$bars.h.mm[k]) else rep(2.5, length(k))
  h[is.na(h)] = 2.5
  has = vapply(k, function(i) { b = edges$bars[[i]]; length(b) > 0 && any(!is.na(b)) }, TRUE)
  ifelse(has, vapply(h + 0.3, .mm2y, 1), 0) ## .mm2y is not vectorized
}

.gt_edge_bars = function(edges, x, y)
{
  if (is.null(edges$bars) || !length(x)) return(invisible(NULL))
  if (!is.list(edges$bars)) edges$bars = as.list(edges$bars)
  done = character(0)
  for (i in seq_along(x))
  {
    b = suppressWarnings(as.numeric(edges$bars[[i]]))
    if (!length(b) || all(is.na(b)) || is.na(x[i]) || is.na(y[i])) next
    key = paste(round(x[i], 2), round(y[i], 4)); if (key %in% done) next; done = c(done, key) ## strand copies of one edge
    w = .mm2x(if (!is.null(edges$bars.w.mm) && !is.na(edges$bars.w.mm[i])) edges$bars.w.mm[i] else 0.7)
    h = .mm2y(if (!is.null(edges$bars.h.mm) && !is.na(edges$bars.h.mm[i])) edges$bars.h.mm[i] else 2.5)
    cl = if (is.null(edges$bars.col)) 'grey40' else { c0 = edges$bars.col[[i]]; if (length(c0) == 1) trimws(strsplit(as.character(c0), ',')[[1]]) else as.character(c0) }
    cl = rep_len(cl, length(b))
    n = length(b); x0 = x[i] - n * w / 2 + (seq_len(n) - 1) * w; y0 = y[i] + .mm2y(0.3)
    graphics::rect(x0, y0, x0 + w * 0.85, y0 + h, col = 'white', border = 'grey80', lwd = 0.3)
    v = pmin(pmax(b, 0), 1); v[is.na(v)] = 0
    graphics::rect(x0, y0, x0 + w * 0.85, y0 + h * v, col = cl, border = NA)
  }
  invisible(NULL)
}
