library(testthat)
library(gTrack)
library(GenomicRanges)

context("raster matrix track")

## synthetic cells x bins copy-number matrix
toy = function(ncell = 300)
{
  bins = GRanges("1", IRanges(seq(1, 99001, 1000), width = 1000))
  m = matrix(2, nrow = length(bins), ncol = ncell, dimnames = list(NULL, paste0("c", 1:ncell)))
  m[40:60, 1:(ncell / 2)] = 8; m[70:80, ] = 1
  bins$M = m
  bins
}

test_that("matrix column and numeric columns both read", {
  b = toy(5)
  expect_equal(dim(gTrack:::.gt_raster_matrix(b, "M")), c(100, 5))
  b2 = granges(b); mcols(b2) = as.data.frame(b$M)
  expect_equal(dim(gTrack:::.gt_raster_matrix(b2, NA)), c(100, 5))
  expect_equal(gTrack:::.gt_raster_colors(c(1, 2, 8), 1, 8, c("blue", "white", "red"), mid = 2), c("#0000FF", "#FFFFFF", "#FF0000"))
})

test_that("raster track draws for thousands of rows, alone and under other tracks; default tracks unaffected", {
  b = toy(2000)
  rt = gTrack(b, raster = TRUE, raster.field = "M", raster.mid = 2, name = "cells", height = 20)
  expect_false(formatting(gTrack(b))$raster)
  segs = GRanges("1", IRanges(c(1, 39001, 61001), c(39000, 61000, 1e5)), y = c(2, 8, 2))
  f = tempfile(fileext = ".pdf"); pdf(f, 5, 4)
  expect_error(suppressWarnings(plot(rt, GRanges("1", IRanges(1, 1e5)))), NA)
  expect_error(suppressWarnings(plot(c(rt, gTrack(segs, y.field = "y", name = "cn")),
                                     c(GRanges("1", IRanges(1, 45000)), GRanges("1", IRanges(55000, 1e5))), style = pub_style())), NA)
  dev.off()
  expect_true(file.size(f) < 2e6) ## one image per window, not one rectangle per cell x bin
})
