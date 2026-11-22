library(testthat)
library(gTrack)
library(GenomicRanges)

context("row bands")

## synthetic toy data only
test_that("band strips, labels and group brackets draw; nothing without band fields", {
  segs = GRanges("1", IRanges(c(1, 2001, 6001), c(2000, 6000, 10000)), y = c(2, 40, 2))
  gt = do.call(c, lapply(1:5, function(k) gTrack(segs, y.field = "y", name = paste0("cell", k))))
  w = GRanges("1", IRanges(1, 10000))
  f = tempfile(fileext = ".pdf"); pdf(f, 4, 4)
  L = suppressWarnings(plot(gt, w, style = pub_style()))
  expect_null(gTrack:::.gt_draw_bands(gt, L))
  formatting(gt)$band.col = c("orange", "orange", "purple", "purple", "purple")
  formatting(gt)$band.label = paste0("cell", 1:5)
  formatting(gt)$band.group = c("CI1", "CI1", "CI2", "CI2", "CI2")
  formatting(gt)$band.group.col = c("orange", "orange", "purple", "purple", "purple")
  expect_error(suppressWarnings(plot(gt, w, style = pub_style(mar.mm = c(6, 25, 2, 2)))), NA)
  expect_error(suppressWarnings(plot(gt, w, band.params = list(offset.mm = 4))), NA)
  dev.off()
})
