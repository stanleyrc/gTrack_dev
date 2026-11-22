library(testthat)
library(gTrack)
library(GenomicRanges)

context("junction support bars")

## synthetic 3-segment graph with one ALT junction carrying per-group support
toy = function(bars = TRUE)
{
  gr = GRanges("1", IRanges(c(1, 3001, 6001), c(3000, 6000, 10000)), strand = "+", y = c(2, 10, 2))
  ed = data.frame(from = c(1, 2, 3), to = c(2, 3, 1), type = c("REF", "REF", "ALT"), not.flat = c(FALSE, FALSE, TRUE),
                  col = c("grey", "grey", "red"), label = c(NA, NA, "J1 (66)"), stringsAsFactors = FALSE)
  if (bars) { ed$bars = list(NULL, NULL, c(1, 0.83, 0.28, 0.17)); ed$bars.col = c(NA, NA, "orange,skyblue,purple,hotpink") }
  gTrack(gr, y.field = "y", edges = ed, name = "graph")
}

test_that("bars draw with classic and bell arcs; tracks without bars unchanged", {
  w = GRanges("1", IRanges(1, 10000))
  f = tempfile(fileext = ".pdf"); pdf(f, 4, 3)
  expect_error(suppressWarnings(plot(toy(), w)), NA)
  expect_error(suppressWarnings(plot(toy(), w, style = pub_style(arc.style = "bell"))), NA)
  expect_error(suppressWarnings(plot(toy(), w, style = pub_style())), NA)
  expect_error(suppressWarnings(plot(toy(FALSE), w, style = pub_style(arc.style = "bell"))), NA)
  dev.off()
  pdf(f, 4, 3); plot.new()
  expect_equal(gTrack:::.gt_bars_h(data.frame(x = 1:2), 1:2), c(0, 0))
  e = data.frame(x = 1:2); e$bars = list(NULL, c(0.5, 1))
  h = gTrack:::.gt_bars_h(e, 1:2); dev.off()
  expect_equal(h[1], 0); expect_true(h[2] > 0)
})
