library(testthat)
library(gTrack)
library(GenomicRanges)

context("per-row y policy")

## synthetic toy data only
toy = function(amp) GRanges("1", IRanges(seq(1, 9901, 100), width = 100), y = c(rep(2, 30), rep(amp, 40), rep(2, 30)))
w = GRanges("1", IRanges(1, 10000))

test_that("gTrack takes y.log1p and y.group; defaults unchanged", {
  g = gTrack(toy(40), y.field = "y", y.log1p = TRUE, y.group = "a")
  expect_true(formatting(g)$y.log1p); expect_equal(formatting(g)$y.group, "a")
  g0 = gTrack(toy(40), y.field = "y")
  expect_false(formatting(g0)$y.log1p); expect_true(is.na(formatting(g0)$y.group))
  o = gTrack:::.gt_y_policy(c(gTrack(toy(40), y.field = "y"), gTrack(toy(400), y.field = "y")), NULL, w)
  expect_true(all(is.na(formatting(o)$y1))) ## own ranges: nothing preset
})

test_that("shared and y.group give one range; user y0/y1 kept; log1p sets the flag", {
  gt = c(gTrack(toy(40), y.field = "y"), gTrack(toy(400), y.field = "y"), gTrack(toy(8), y.field = "y", y1 = 10))
  o = gTrack:::.gt_y_policy(gt, "shared", w)
  expect_equal(formatting(o)$y1[1], formatting(o)$y1[2])
  expect_true(formatting(o)$y1[1] >= 400)
  expect_equal(formatting(o)$y1[3], 10)
  gt2 = c(gTrack(toy(40), y.field = "y", y.group = "g"), gTrack(toy(400), y.field = "y", y.group = "g"), gTrack(toy(8), y.field = "y"))
  o2 = gTrack:::.gt_y_policy(gt2, NULL, w)
  expect_equal(formatting(o2)$y1[1], formatting(o2)$y1[2]); expect_true(is.na(formatting(o2)$y1[3]))
  expect_true(all(formatting(gTrack:::.gt_y_policy(gt, "log1p", w))$y.log1p))
  expect_error(gTrack:::.gt_y_policy(gt, "bogus", w))
  f = tempfile(fileext = ".pdf"); pdf(f, 4, 3)
  for (p in list(NULL, "own", "shared", "log1p", "shared.log1p"))
    expect_error(suppressWarnings(plot(gt, w, y.policy = p, style = pub_style())), NA)
  dev.off()
})
