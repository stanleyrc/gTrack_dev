## Per-track loose-end toggle (additive; default draws loose ends as before).
##
## gGraph tracks (gGnome gg$gtrack()) hold loose ends as extra ranges with $loose == TRUE, joined to their node by
## edges of type "loose". Setting formatting(gt)$loose = FALSE for a track (or plot(gt, loose = FALSE) for all tracks)
## removes the edges to those ranges before drawing (the blue loose-end marks), e.g. the loose-end dots where a consensus graph is cut at the
## window edge.
.gt_drop_loose = function(.Object, loose.arg = NULL)
{
  ## the loose ranges themselves are kept (gGnome draws them transparent, and they take part in the automatic segment
  ## thickness), only the edges to them are removed, so the rest of the track looks exactly as before
  f = formatting(.Object)
  lz = if (!is.null(loose.arg)) rep_len(as.logical(loose.arg), length(.Object))
       else if (!is.null(f$loose)) as.logical(f$loose) else rep(NA, length(.Object))
  for (j in which(lz %in% FALSE))
  {
    ed = .Object@edges[[j]]
    if (is.null(ed) || NROW(ed) == 0) next
    ed = as.data.frame(ed, stringsAsFactors = FALSE)
    drop = if (!is.null(ed$type)) ed$type %in% 'loose' else rep(FALSE, nrow(ed))
    d = .Object@data[[j]]
    if (inherits(d, 'GRanges') && !is.null(d$loose) && all(c('from', 'to') %in% names(ed)))
      drop = drop | d$loose[ed$from] %in% TRUE | d$loose[ed$to] %in% TRUE
    if (any(drop)) .Object@edges[[j]] = ed[!drop, , drop = FALSE]
  }
  .Object
}
