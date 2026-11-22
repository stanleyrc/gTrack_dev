library(testthat)
library(gTrack)
library(GenomicRanges)

context("pub_style title.pos")

## synthetic toy data only
test_that("title.pos top (default), left and inside all draw", {
  expect_equal(pub_style()$title.pos, "top")
  segs = GRanges("1", IRanges(c(1, 2001, 6001), c(2000, 6000, 10000)), y = c(2, 200, 2))
  gt = c(gTrack(segs, y.field = "y", name = "cov", height = 4), gTrack(segs, y.field = "y", name = "graph", height = 4))
  f = tempfile(fileext = ".pdf"); pdf(f, 4, 3)
  for (tp in c("top", "left", "inside"))
    expect_error(suppressWarnings(plot(gt, GRanges("1", IRanges(1, 10000)), style = pub_style(title.pos = tp))), NA)
  formatting(gt)$name.pos = c("left", "inside") ## per-track override
  expect_error(suppressWarnings(plot(gt, GRanges("1", IRanges(1, 10000)), style = pub_style())), NA)
  dev.off()
})
