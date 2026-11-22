## Side-by-side triangles on one color scale, and diverging difference maps (additive).

.gt_tri_shared = function(.Object, tri.shared, windows)
{
  f = formatting(.Object)
  js = which(f$triangle %in% TRUE)
  if (!is.null(f$tri.bg.value)) js = js[is.na(f$tri.bg.value[js])] ## difference maps (gt_tri_diff) keep their own diverging scale
  if (length(js) < 1) return(.Object)
  if (is.numeric(tri.shared) && length(tri.shared) == 2) lim = sort(tri.shared)
  else if (all(!is.na(f$cmap.min[js])) && all(!is.na(f$cmap.max[js]))) lim = c(min(f$cmap.min[js]), max(f$cmap.max[js]))
  else
  {
    v = unlist(lapply(js, function(j) {
      m = .Object@mdata[[j]]; d = .Object@data[[j]]
      if (!inherits(d, 'GRanges') || is.null(dim(m)) || nrow(m) != length(d)) return(NULL)
      ix = which(gUtils::gr.in(d, windows))
      if (!length(ix)) return(NULL)
      mm = m[ix, ix, drop = FALSE]
      if (inherits(mm, 'Matrix')) mm = as(mm, 'CsparseMatrix')@x else mm = as.vector(mm)
      mm[is.finite(mm) & mm != 0]
    }))
    lim = if (length(v)) as.numeric(stats::quantile(v, c(0.01, 0.99))) else c(0, 1)
    if (any(!is.na(f$cmap.min[js]))) lim[1] = min(c(lim[1], f$cmap.min[js]), na.rm = TRUE)
    if (any(!is.na(f$cmap.max[js]))) lim[2] = max(c(lim[2], f$cmap.max[js]), na.rm = TRUE)
  }
  formatting(.Object)$cmap.min[js] = lim[1]
  formatting(.Object)$cmap.max[js] = lim[2]
  formatting(.Object)$legend[js] = FALSE
  formatting(.Object)$legend[max(js)] = TRUE ## one color bar, on the top triangle track
  .Object
}

#' @name gt_tri_diff
#' @title Difference map of two triangle (mdata) tracks
#' @description
#' Returns a triangle gTrack of mdata(a) - mdata(b) on a diverging palette centered on 0 (0 is the white background and is
#' not drawn per pixel). Both tracks must be built on the same bins (e.g. two Hi-C clusters from one gMatrix binning).
#' @param a,b length-1 gTracks with mdata on identical ranges
#' @param name track name
#' @param lim color limit; the scale is c(-lim, lim). Default: 99th percentile of |difference|
#' @param cols diverging palette (low, mid, high)
#' @param ... further gTrack formatting
#' @return gTrack
#' @export
gt_tri_diff = function(a, b, name = 'difference', lim = NA, cols = c('#2166AC', 'white', '#B2182B'), ...)
{
  ga = a@data[[1]]; gb = b@data[[1]]
  if (!inherits(ga, 'GRanges') || length(ga) != length(gb) || !all(ga == gb))
    stop('gt_tri_diff: both tracks need mdata on identical ranges')
  d = a@mdata[[1]] - b@mdata[[1]]
  if (inherits(d, 'Matrix')) d = Matrix::drop0(d)
  if (is.na(lim))
  {
    x = if (inherits(d, 'Matrix')) as(d, 'CsparseMatrix')@x else as.vector(d)
    lim = if (length(x)) as.numeric(stats::quantile(abs(x[is.finite(x)]), 0.99)) else 1
    if (!is.finite(lim) || lim <= 0) lim = 1
  }
  gt = gTrack(granges(ga), mdata = d, name = name, cmap.min = -lim, cmap.max = lim, ...)
  colormap(gt) = list(value = cols)
  formatting(gt)$tri.bg.value = 0
  gt
}
