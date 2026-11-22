library(testthat)
library(gTrack)
library(GenomicRanges)

context("loose-end toggle")

## synthetic gGraph-like track: 3 nodes + 2 loose ends (as gGnome gg$gtrack() builds them)
toy_graph = function()
{
  gr = GRanges("1", IRanges(c(1, 3001, 6001, 1, 10000), c(3000, 6000, 10000, 1, 10000)), strand = "+",
               y = c(2, 6, 2, 2, 2), loose = c(FALSE, FALSE, FALSE, TRUE, TRUE))
  ed = data.frame(from = c(1, 2, 4, 3), to = c(2, 3, 1, 5), type = c("REF", "REF", "loose", "loose"),
                  col = c("grey", "grey", "blue", "blue"), stringsAsFactors = FALSE)
  gTrack(gr, y.field = "y", edges = ed, name = "graph")
}

test_that("loose = FALSE drops loose ranges and edges, default keeps them", {
  gt = toy_graph()
  kept = gTrack:::.gt_drop_loose(gt)
  expect_equal(length(kept@data[[1]]), 5)
  formatting(gt)$loose = FALSE
  dropped = gTrack:::.gt_drop_loose(gt)
  expect_equal(length(dropped@data[[1]]), 5) ## ranges kept (transparent in gGnome), edges dropped
  expect_equal(nrow(dropped@edges[[1]]), 2)
  expect_true(all(dropped@edges[[1]]$type == "REF"))
  expect_equal(nrow(gTrack:::.gt_drop_loose(toy_graph(), FALSE)@edges[[1]]), 2) ## plot(gt, loose = FALSE)
  f = tempfile(fileext = ".pdf"); pdf(f, 4, 3)
  expect_error(suppressWarnings(plot(toy_graph(), GRanges("1", IRanges(1, 10000)), loose = FALSE)), NA)
  expect_error(suppressWarnings(plot(gt, GRanges("1", IRanges(1, 10000)))), NA)
  dev.off()
})
