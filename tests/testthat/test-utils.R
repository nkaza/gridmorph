test_that(".valid_cells excludes NA and 0, keeps everything else", {
    rast <- terra::rast(nrows = 1, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 1, crs = "local")
    terra::values(rast) <- c(1, 5, 0, NA)
    expect_equal(as.logical(terra::values(.valid_cells(rast))), c(TRUE, TRUE, FALSE, FALSE))
})

test_that(".valid_cells errors on non-SpatRaster input", {
    expect_error(.valid_cells(1:5), "SpatRaster")
})

test_that(".mass_raster returns rast itself when weighted, a constant when not", {
    rast <- terra::rast(nrows = 1, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 1, crs = "local")
    terra::values(rast) <- c(3, 7, 0, NA)
    valid <- .valid_cells(rast)

    w_true <- .mass_raster(rast, valid, weighted = TRUE)
    expect_equal(as.vector(terra::values(w_true)), c(3, 7, 0, NA))

    w_false <- .mass_raster(rast, valid, weighted = FALSE)
    expect_equal(as.vector(terra::values(w_false))[c(1, 2)], c(1, 1))
    expect_true(all(is.na(as.vector(terra::values(w_false))[c(3, 4)])))
})

test_that(".mass_raster warns on a categorical raster with weighted = TRUE, not with weighted = FALSE", {
    rast <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1, crs = "local")
    terra::values(rast) <- c(1, 3, 0)
    levels(rast) <- data.frame(id = c(0, 1, 3), category = c("bg", "classA", "classB"))
    valid <- .valid_cells(rast)

    expect_warning(.mass_raster(rast, valid, weighted = TRUE), "categorical")
    expect_no_warning(.mass_raster(rast, valid, weighted = FALSE))
})

test_that(".mass_raster does not warn on a plain numeric raster, even with few distinct integer values", {
    rast <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1, crs = "local")
    terra::values(rast) <- c(1, 3, 0)  # same values as above, but NOT marked categorical via levels<-
    valid <- .valid_cells(rast)
    expect_no_warning(.mass_raster(rast, valid, weighted = TRUE))
})

test_that(".sample_valid_points draws WITH replacement (a genuine i.i.d. requirement)", {
    r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "local")
    terra::values(r) <- 0
    r[2:4, 2:4] <- 1
    valid <- .valid_cells(r)
    w <- terra::ifel(valid, 1, NA)
    # 9 valid cells; requesting far more than that only works if sampling
    # is WITH replacement, not without - the defining i.i.d. requirement
    pts <- .sample_valid_points(valid, w, size = 500)
    expect_equal(nrow(pts), 500)
    expect_lt(nrow(unique(pts)), 500)  # necessarily some repeats
    expect_equal(nrow(unique(pts)), 9)  # every valid cell gets drawn at least once
})

test_that(".sample_valid_points's empirical frequencies match the target weight distribution", {
    r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "local")
    terra::values(r) <- 0
    r[2:4, 2:4] <- 1
    valid <- .valid_cells(r)
    w <- r
    wv <- rep(NA_real_, terra::ncell(r))
    wv[as.vector(terra::values(valid))] <- c(10, 1, 1, 1, 1, 1, 1, 1, 1)  # one cell 10x the rest
    terra::values(w) <- wv

    set.seed(1)
    pts <- .sample_valid_points(valid, w, size = 50000)
    key <- paste(pts[, 1], pts[, 2])
    prop <- table(key) / nrow(pts)
    # the heavy cell should sit at ~10/18 = 0.5556; every other cell at ~1/18 = 0.0556
    expect_equal(unname(max(prop)), 10 / 18, tolerance = 0.02)
    expect_equal(unname(min(prop)), 1 / 18, tolerance = 0.02)
})

test_that(".sample_valid_points's draw sequence shows no lag-1 autocorrelation (independence, not just matching marginals)", {
    r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "local")
    terra::values(r) <- 0
    r[2:4, 2:4] <- 1
    valid <- .valid_cells(r)
    w <- terra::ifel(valid, 1, NA)

    set.seed(1)
    pts <- .sample_valid_points(valid, w, size = 20000)
    idx <- as.integer(factor(paste(pts[, 1], pts[, 2])))
    ac <- cor(idx[-length(idx)], idx[-1])
    expect_lt(abs(ac), 0.03)
})

make_large_gradient_raster <- function(n = 100) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    terra::values(r) <- 1
    list(r = r, valid = .valid_cells(r), w = terra::init(r, "x") + terra::init(r, "y") + 1)
}

test_that(".sample_valid_points takes the binned path once n_valid exceeds n_bins, not before", {
    g <- make_large_gradient_raster()  # 10000 valid cells > default n_bins = 1000
    b <- .bin_index_raster(g$valid, g$w, 1000)
    expect_equal(b$mode, "binned")
    pts <- .sample_valid_points(g$valid, g$w, size = 2000)
    expect_equal(nrow(pts), 2000)
})

test_that("the binned path's empirical mean matches the density-weighted target, not the plain mean", {
    g <- make_large_gradient_raster()
    wvals <- as.vector(terra::values(g$w))
    target_mean <- sum(wvals^2) / sum(wvals)  # E[W] under density proportional to W

    set.seed(1)
    pts <- .sample_valid_points(g$valid, g$w, size = 50000)
    sampled_w <- terra::extract(g$w, pts)[[1]]
    expect_equal(mean(sampled_w), target_mean, tolerance = 0.02)
    expect_false(isTRUE(all.equal(mean(sampled_w), mean(wvals), tolerance = 0.02)))
})

test_that("the binned path's draw sequence shows no lag-1 autocorrelation (the post-bin shuffle works)", {
    g <- make_large_gradient_raster()
    set.seed(1)
    pts <- .sample_valid_points(g$valid, g$w, size = 20000)
    idx <- as.integer(factor(paste(pts[, 1], pts[, 2])))
    ac <- cor(idx[-length(idx)], idx[-1])
    expect_lt(abs(ac), 0.03)
})

test_that("a bin with very few member cells doesn't crash the binned sampling path", {
    # regression test for a real terra bug hit during development:
    # terra::spatSample(method = "random") on a mask with exactly one
    # matching cell crashes inside terra:::.sampleCellsMemory() ("length(n)
    # == 1L is not TRUE") on a small/single-block raster - .sample_from_bins()
    # was rewritten to draw cell indices and read values directly instead
    # of delegating to spatSample(method = "random") for exactly this reason
    n <- 40
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    terra::values(r) <- 1
    valid <- .valid_cells(r)
    # a near-constant weight with ONE deliberately isolated outlier cell,
    # forcing a sparsely-populated top bin once n_bins forces binning
    wv <- rep(1, terra::ncell(r))
    wv[1] <- 1000
    w <- r
    terra::values(w) <- wv

    expect_no_error(pts <- .sample_valid_points(valid, w, size = 500, n_bins = 100))
    expect_equal(nrow(pts), 500)
})

test_that("degenerate (constant-weight) mode routes through method = \"random\", memory-safe regardless of n_valid_cells", {
    # regression test for a real gap found during development: an earlier
    # version routed "degenerate" mode through method = "weights" alongside
    # "exact" mode, which is only cheap when n_valid is ALSO small - for a
    # constant weight on a large raster (e.g. an ordinary big binary mask),
    # that reintroduces the exact n_valid-scaling as.data.frame() cost this
    # whole two-stage design exists to avoid. This doesn't test memory
    # directly, but confirms the degenerate branch is taken and produces
    # correct output regardless of raster size relative to n_bins.
    g <- make_large_gradient_raster()  # 10000 cells, > default n_bins
    w_const <- terra::ifel(g$valid, 1, NA)
    b <- .bin_index_raster(g$valid, w_const, 1000)
    expect_equal(b$mode, "degenerate")
    pts <- .sample_valid_points(g$valid, w_const, size = 2000)
    expect_equal(nrow(pts), 2000)
})

#' gm_convexity_index()/gm_span_index() split .sample_valid_points()'s output
#' into consecutive (x1, x2) pairs - LINES, not just points, need to be
#' i.i.d.: each pair's own two endpoints independent of each other, AND
#' different pairs independent of each other. Points being i.i.d. is
#' necessary but not sufficient for this on its own if anything downstream
#' groups or reorders them before pairing - which is exactly what
#' .sample_from_bins() does (assembles results one bin at a time) before
#' its own explicit shuffle. This checks the PAIRING itself, using a
#' spatially-correlated weight (a left-to-right gradient) specifically
#' because that's the adversarial case: if pairing were NOT properly
#' randomised, paired points would spuriously cluster together in space
#' (same-bin cells tend to be spatially near each other for a smooth
#' gradient), which a plain "points are i.i.d." check would never catch.
test_that("consecutive pairs formed from the sample show no spurious within-pair correlation", {
    check_pair_independence <- function(valid, w, size, n_bins) {
        pts <- .sample_valid_points(valid, w, size, n_bins = n_bins)
        n_pairs <- nrow(pts) %/% 2
        x1 <- pts[seq(1, 2 * n_pairs, 2), , drop = FALSE]
        x2 <- pts[seq(2, 2 * n_pairs, 2), , drop = FALSE]
        # under true independence, corr(x1_x, x2_x) is ~0 REGARDLESS of how
        # strongly x correlates with weight - a real pairing bug would show
        # up as a spurious positive correlation here, not as "close to 0"
        cor(x1[, 1], x2[, 1])
    }

    r <- terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 100, ymin = 0, ymax = 100, crs = "local")
    terra::values(r) <- 1
    valid <- .valid_cells(r)
    w <- terra::init(r, "x")  # a strong left-to-right gradient

    set.seed(1)
    # n_valid = 10000 here; n_bins = 20 forces heavy binning (coarse bins,
    # the scenario most likely to expose a broken shuffle) alongside the
    # default n_bins = 1000 (light binning) for the same raster
    expect_lt(abs(check_pair_independence(valid, w, size = 20000, n_bins = 1000)), 0.03)
    expect_lt(abs(check_pair_independence(valid, w, size = 20000, n_bins = 20)), 0.03)
})
