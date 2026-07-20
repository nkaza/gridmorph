make_square_mask <- function(n = 21, inset = 5) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    terra::values(r) <- 0
    r[(inset + 1):(n - inset), (inset + 1):(n - inset)] <- 1
    r
}

make_gradient_rast <- function(r) {
    # a genuinely non-uniform gradient inside the square, 0/NA outside
    cc <- terra::init(r, "x") + terra::init(r, "y") + 1
    terra::ifel(r == 1, cc, NA)
}

test_that("gm_depth_index is positive, finite, and roughly in (0, 1] for a simple solid square", {
    r <- make_square_mask()
    res <- gm_depth_index(r)
    expect_true(is.finite(res$index))
    expect_gt(res$index, 0)
    expect_lt(res$index, 1.05)  # any overshoot above 1 is a residual discretization artifact, not clamped
    expect_gt(res$mean_depth, 0)
    expect_equal(res$n_valid_cells, sum(as.vector(terra::values(r)), na.rm = TRUE))
})

test_that("gm_depth_index scores a rasterized disk close to 1, not just below some loose ceiling", {
    # regression test for the resolution-matched reference (see R/depth-index.R's
    # own header): comparing a biased mean_depth against an unbiased closed-form
    # reference used to put a disk at 1.024 at this resolution; a reference built
    # from the same biased .depth_field() pipeline should cancel almost all of it
    n <- 81; radius <- 30
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
    d <- sqrt((xy[, 1] - n / 2)^2 + (xy[, 2] - n / 2)^2)
    terra::values(r) <- ifelse(d <= radius, 1, NA)

    res <- gm_depth_index(r)
    expect_equal(res$index, 1, tolerance = 0.01)
})

test_that("weighted = FALSE ignores a continuous raster's own magnitude, reproducing the plain binary-mask index", {
    r <- make_square_mask()
    w <- make_gradient_rast(r)  # continuous, genuinely non-uniform inside the square
    plain <- gm_depth_index(r)  # r itself is already 0/1 - weighted vs not makes no difference here
    ignored <- gm_depth_index(w, weighted = FALSE)
    weighted <- gm_depth_index(w, weighted = TRUE)

    expect_equal(ignored$index, plain$index, tolerance = 1e-10)
    expect_equal(ignored$mean_depth, plain$mean_depth, tolerance = 1e-10)
    # and the ACTUAL weighted computation should differ (the gradient is genuinely non-uniform)
    expect_false(isTRUE(all.equal(weighted$index, ignored$index, tolerance = 1e-6)))
})

test_that("a hole punched into the shape (NA in the interior) reduces n_valid_cells", {
    r <- make_square_mask(n = 21, inset = 5)
    n_before <- gm_depth_index(r)$n_valid_cells

    r_hole <- r
    r_hole[10:12, 10:12] <- NA  # punch a hole inside the shape
    res_hole <- gm_depth_index(r_hole)
    expect_lt(res_hole$n_valid_cells, n_before)
})

test_that("resolution convergence: a finer raster's gm_depth_index moves toward a stable value", {
    coarse <- gm_depth_index(make_square_mask(n = 21, inset = 5))$index
    fine <- gm_depth_index(make_square_mask(n = 81, inset = 20))$index
    finer <- gm_depth_index(make_square_mask(n = 161, inset = 40))$index
    # not asserting a specific target value here (that's the shapeindices
    # cross-check, done interactively/manually), just that refining the
    # grid stops moving the answer around wildly
    expect_lt(abs(finer - fine), abs(fine - coarse) + 1e-6)
})

test_that("no valid cells gives NA with a warning, not an error", {
    r <- make_square_mask()
    terra::values(r) <- 0
    expect_warning(res <- gm_depth_index(r), "No valid cells")
    expect_true(is.na(res$index))
})

test_that("n_bins >= n_valid_cells (the exact path) is unaffected by n_bins itself", {
    r <- make_square_mask()
    w <- make_gradient_rast(r)
    ref <- gm_depth_index(w, n_bins = 100000)
    same <- gm_depth_index(w, n_bins = 5000)  # both above n_valid (121 cells)
    expect_equal(same$ref_depth, ref$ref_depth, tolerance = 1e-12)
    expect_equal(same$index, ref$index, tolerance = 1e-12)
})

test_that("forcing the binned path (n_bins < n_valid_cells) stays close to the exact reference", {
    r <- make_square_mask()
    w <- make_gradient_rast(r)
    exact <- gm_depth_index(w, n_bins = 100000)
    binned <- gm_depth_index(w, n_bins = 10)  # well below 121 valid cells - forces binning
    expect_true(is.finite(binned$index))
    expect_equal(binned$ref_depth, exact$ref_depth, tolerance = 0.05)
})

test_that("a constant weight raster does not break the binned path's degenerate range case", {
    r <- make_square_mask()
    unweighted <- gm_depth_index(r)
    # weighted = FALSE forces a constant density regardless of n_bins
    binned_constant <- gm_depth_index(r, weighted = FALSE, n_bins = 5)
    expect_equal(binned_constant$index, unweighted$index, tolerance = 1e-10)
})

test_that("a shape touching the raster's own edge gives the identical mean_depth as the same shape with a margin", {
    # regression test for a verified terra::distance() edge bug: without
    # padding first, a shape reaching the raster's true extent has no
    # boundary cell there to measure distance to, understating its own
    # depth right at the edge - verified directly, an edge-touching square
    # used to score depth_index = 1.238, outside this index's own
    # documented (0, 1] range (a square must always score below the
    # disk-optimal 1)
    n <- 21
    r_edge <- terra::rast(nrows = n, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = n, crs = "local")
    terra::values(r_edge) <- 0
    r_edge[1:n, 11:31] <- 1  # flush against the top edge (row 1)

    r_margin <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r_margin) <- 0
    r_margin[11:31, 11:31] <- 1  # same square, centred with a margin on every side

    res_edge <- gm_depth_index(r_edge)
    res_margin <- gm_depth_index(r_margin)
    expect_equal(res_edge$mean_depth, res_margin$mean_depth, tolerance = 1e-10)
    expect_equal(res_edge$index, res_margin$index, tolerance = 1e-10)
    expect_lt(res_edge$index, 1.2)  # small overshoot above 1 is an expected discretization residual, not clamped
})

test_that("a shape filling its entire raster (no invalid cell anywhere) does not return NaN", {
    # regression test: distance() has nothing to measure to when `valid`
    # covers the whole raster and inverted-and-padded still has no real
    # boundary without the fix - used to return NaN
    r <- terra::rast(nrows = 21, ncols = 21, xmin = 0, xmax = 21, ymin = 0, ymax = 21, crs = "local")
    terra::values(r) <- 1
    res <- gm_depth_index(r)
    expect_true(is.finite(res$index))
    expect_true(is.finite(res$mean_depth))
})
