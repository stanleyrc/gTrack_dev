## Per-row y policy (additive; the default y.policy = NULL keeps every track on its own range as before).
.gt_y_policy = function(.Object, y.policy = NULL, windows)
{
  if (!is.null(y.policy) && !y.policy %in% c('own', 'shared', 'log1p', 'shared.log1p'))
    stop("y.policy must be NULL, 'own', 'shared', 'log1p' or 'shared.log1p'")
  f = formatting(.Object)
  hasy = !is.na(f$y.field)
  if (!any(hasy)) return(.Object)
  if (is.null(f$y.group)) formatting(.Object)$y.group = NA
  if (is.null(f$y.log1p)) formatting(.Object)$y.log1p = FALSE
  formatting(.Object)$y.log1p[is.na(formatting(.Object)$y.log1p)] = FALSE
  if (!is.null(y.policy) && grepl('log1p', y.policy))
    formatting(.Object)$y.log1p[hasy] = TRUE
  if (!is.null(y.policy) && grepl('shared', y.policy))
  {
    grp = as.character(formatting(.Object)$y.group)
    grp[hasy & is.na(grp)] = '.shared'
    formatting(.Object)$y.group = grp
  }
  grp = as.character(formatting(.Object)$y.group)
  for (g in setdiff(unique(grp[hasy]), NA))
  {
    js = which(hasy & grp %in% g)
    if (length(js) < 2) next
    lim = sapply(js, function(j) {
      tt = tryCatch(extract_data_from_tmp_dat(.Object, j, windows), error = function(e) NULL)
      if (is.null(tt) || !inherits(tt$t, c('GRanges', 'GRangesList')) || !(f$y.field[j] %in% names(values(tt$t)))) return(c(NA, NA))
      o = tryCatch(format_yfield_limits(tt$o, j, tt$t, FALSE, tt$w), error = function(e) NULL)
      if (is.null(o)) c(NA, NA) else c(formatting(o)$y0[j], formatting(o)$y1[j])
    })
    if (all(is.na(lim))) next
    y0 = min(lim[1, ], na.rm = TRUE); y1 = max(lim[2, ], na.rm = TRUE)
    for (j in js)
    {
      if (is.na(formatting(.Object)$y0[j])) formatting(.Object)$y0[j] = y0
      if (is.na(formatting(.Object)$y1[j])) formatting(.Object)$y1[j] = y1
    }
  }
  .Object
}
