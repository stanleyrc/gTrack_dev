library(testthat)
library(gTrack)
library(GenomicRanges)

context("layout return value")

## synthetic toy data only
test_that("plot.gTrack returns the track and window layout invisibly", {
  segs = GRanges("1", IRanges(c(1, 2001, 6001), c(2000, 6000, 10000)), y = c(2, 40, 2))
  gt = c(gTrack(segs, y.field = "y", name = "A"), gTrack(segs, y.field = "y", name = "B", height = 20))
  f = tempfile(fileext = ".pdf"); pdf(f, 4, 3)
  vis = withVisible(plot(gt, c(GRanges("1", IRanges(1, 4000)), GRanges("1", IRanges(5000, 10000)))))
  L = vis$value
  dev.off()
  expect_false(vis$visible)
  expect_s3_class(L, "gTrackLayout")
  expect_equal(nrow(L$tracks), 2)
  expect_equal(L$tracks$name, c("A", "B"))
  expect_true(L$tracks$y1[2] > L$tracks$y1[1]) ## track 2 drawn above track 1
  expect_true(all(L$tracks$ndc.y0 >= 0 & L$tracks$ndc.y1 <= 1))
  expect_equal(nrow(L$windows), 2)
  x = gt_layout_x(L, c("1", "1", "1"), c(1, 4000, 4500))
  expect_true(x[1] < x[2]); expect_true(is.na(x[3])) ## 4500 is in the gap between windows
  y = gt_layout_y(L, "A", c(L$tracks$ydata0[1], L$tracks$ydata1[1]))
  expect_equal(y, c(L$tracks$y0[1], L$tracks$y1[1]))
})
