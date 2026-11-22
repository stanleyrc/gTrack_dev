library(testthat)
library(gTrack)
library(GenomicRanges)

context("shared triangle color scale and difference maps")

## synthetic contact matrices on 50 bins: B has twice the contacts of A plus an extra block
toy = function()
{
  bins = GRanges("1", IRanges(seq(1, 49001, 1000), width = 1000))
  d = abs(outer(1:50, 1:50, "-")); A = 100 / (1 + d); B = 2 * A; B[20:30, 20:30] = B[20:30, 20:30] + 80
  list(bins = bins, A = A, B = B)
}

test_that("tri.shared gives one range and one legend; default unchanged", {
  x = toy(); w = GRanges("1", IRanges(1, 50000))
  gt = c(gTrack(x$bins, mdata = x$A, name = "A", cmap.min = 0, cmap.max = 50), gTrack(x$bins, mdata = x$B, name = "B", cmap.min = 0, cmap.max = 150))
  o = gTrack:::.gt_tri_shared(gt, TRUE, w)
  expect_equal(formatting(o)$cmap.max, c(150, 150)); expect_equal(formatting(o)$legend, c(FALSE, TRUE))
  expect_equal(formatting(gTrack:::.gt_tri_shared(gt, c(0, 80), w))$cmap.max, c(80, 80))
  gt2 = c(gTrack(x$bins, mdata = x$A, name = "A"), gTrack(x$bins, mdata = x$B, name = "B"))
  o2 = gTrack:::.gt_tri_shared(gt2, TRUE, w)
  expect_equal(formatting(o2)$cmap.max[1], formatting(o2)$cmap.max[2])
  f = tempfile(fileext = ".pdf"); pdf(f, 4, 5)
  expect_error(suppressWarnings(plot(gt, w)), NA)
  expect_error(suppressWarnings(plot(gt, w, tri.shared = TRUE)), NA)
  dev.off()
})

test_that("gt_tri_diff builds a diverging difference track", {
  x = toy(); w = GRanges("1", IRanges(1, 50000))
  a = gTrack(x$bins, mdata = x$A); b = gTrack(x$bins, mdata = x$B)
  d = gt_tri_diff(b, a, name = "B - A")
  expect_equal(formatting(d)$cmap.min, -formatting(d)$cmap.max)
  expect_equal(formatting(d)$tri.bg.value, 0)
  expect_equal(unname(mdata(d)[[1]][25, 25]), x$B[25, 25] - x$A[25, 25])
  expect_error(gt_tri_diff(a, gTrack(x$bins[1:10], mdata = x$A[1:10, 1:10])))
  o = gTrack:::.gt_tri_shared(c(d, a, b), TRUE, w) ## the difference track keeps its own scale
  expect_equal(formatting(o)$cmap.min[1], formatting(d)$cmap.min)
  f = tempfile(fileext = ".pdf"); pdf(f, 4, 3)
  expect_error(suppressWarnings(plot(d, w)), NA)
  dev.off()
})
