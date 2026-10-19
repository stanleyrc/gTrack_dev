library(testthat)
library(gTrack)
library(GenomicRanges)

context("pub_style")

## synthetic toy data only
toy = function()
{
  set.seed(1)
  segs = GRanges("1", IRanges(c(1, 2001, 6001), c(2000, 6000, 10000)), y = c(2, 40, 2))
  pts = GRanges("1", IRanges(seq(1, 9901, 100), width = 1), y = pmax(0, c(rep(2, 20), rep(40, 40), rep(2, 40)) + rnorm(100)))
  list(segs = segs, pts = pts)
}

test_that("pub_style returns a style list", {
  st = pub_style(fontsize = 6)
  expect_true(inherits(st, "gTrackStyle"))
  expect_equal(st$fontsize, 6)
  expect_equal(length(st$mar.mm), 4)
})

test_that("plot.gTrack draws with and without a style, one or several windows", {
  x = toy()
  gt = c(gTrack(x$pts, y.field = "y", circles = TRUE, name = "coverage"), gTrack(x$segs, y.field = "y", name = "segments"))
  f = tempfile(fileext = ".pdf"); pdf(f, 4, 3)
  expect_silent(suppressWarnings(plot(gt, GRanges("1", IRanges(1, 10000)), style = pub_style())))
  expect_silent(suppressWarnings(plot(gt, c(GRanges("1", IRanges(1, 4000)), GRanges("1", IRanges(5000, 10000))), style = pub_style())))
  suppressWarnings(plot(gt, GRanges("1", IRanges(1, 10000))))
  dev.off()
  expect_true(file.size(f) > 0)
  expect_null(getOption("gTrack.style")) ## style is not left behind
})

test_that("gtrack_grob returns a grob", {
  skip_if_not_installed("gridGraphics")
  x = toy()
  g = gtrack_grob(gTrack(x$segs, y.field = "y"), GRanges("1", IRanges(1, 10000)))
  expect_true(inherits(g, "grob"))
})
