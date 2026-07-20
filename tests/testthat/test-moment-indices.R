make_square_mask <- function(n = 41, inset = 15) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    terra::values(r) <- 0
    r[(inset + 1):(n - inset), (inset + 1):(n - inset)] <- 1
    r
}

make_gradient_rast <- function(r) {
    cc <- terra::init(r, "x") + terra::init(r, "y") + 1  # genuinely non-uniform, 0/NA outside
    terra::ifel(r == 1, cc, NA)
}

test_that("gm_moment_of_inertia_index is in (0, 1] for a simple square, high but not 1", {
    r <- make_square_mask()
    res <- gm_moment_of_inertia_index(r)
    expect_true(is.finite(res$index))
    expect_gt(res$index, 0)
    expect_lt(res$index, 1.05)
})

test_that("gm_moment_isotropy_index is close to 1 for a square (4-fold symmetric)", {
    r <- make_square_mask()
    res <- gm_moment_isotropy_index(r)
    expect_equal(res$index, 1, tolerance = 0.02)
})

test_that("gm_directional_balance_index is close to 1 for a centered symmetric square", {
    r <- make_square_mask()
    res <- gm_directional_balance_index(r)
    expect_equal(res$index, 1, tolerance = 0.02)
    expect_equal(res$R, 0, tolerance = 0.02)
})

test_that("gm_moment_isotropy_index is low for an elongated rectangle", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[18:23, 2:39] <- 1  # long thin horizontal rectangle
    res <- gm_moment_isotropy_index(r)
    expect_lt(res$index, 0.2)
})

test_that("weighted = FALSE ignores a continuous raster's own magnitude, all three indices", {
    r <- make_square_mask()
    w <- make_gradient_rast(r)

    moi_plain <- gm_moment_of_inertia_index(r); moi_ignored <- gm_moment_of_inertia_index(w, weighted = FALSE)
    iso_plain <- gm_moment_isotropy_index(r); iso_ignored <- gm_moment_isotropy_index(w, weighted = FALSE)
    db_plain <- gm_directional_balance_index(r); db_ignored <- gm_directional_balance_index(w, weighted = FALSE)

    expect_equal(moi_ignored$index, moi_plain$index, tolerance = 1e-6)
    expect_equal(iso_ignored$index, iso_plain$index, tolerance = 1e-6)
    expect_equal(db_ignored$index, db_plain$index, tolerance = 1e-6)

    # and the ACTUAL weighted computation differs (the gradient is genuinely non-uniform)
    moi_weighted <- gm_moment_of_inertia_index(w, weighted = TRUE)
    expect_false(isTRUE(all.equal(moi_weighted$index, moi_ignored$index, tolerance = 1e-6)))
})

test_that("multi-part shape is treated as one combined shape, not per-patch", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[5:10, 5:10] <- 1
    r[30:35, 30:35] <- 1   # a second, disjoint patch
    res <- gm_moment_of_inertia_index(r)
    expect_true(is.finite(res$index))
    # the combined centroid should sit between the two patches, not inside either
    expect_gt(res$centroid["x"], 10)
    expect_lt(res$centroid["x"], 30)
})

test_that("no valid cells gives NA with a warning for all three indices", {
    r <- make_square_mask(); terra::values(r) <- 0
    expect_warning(res1 <- gm_moment_of_inertia_index(r), "No valid cells")
    expect_warning(res2 <- gm_moment_isotropy_index(r), "No valid cells")
    expect_warning(res3 <- gm_directional_balance_index(r), "No valid cells")
    expect_true(is.na(res1$index)); expect_true(is.na(res2$index)); expect_true(is.na(res3$index))
})

test_that("n_bins >= n_valid_cells (the exact path) is unaffected by n_bins itself", {
    r <- make_square_mask()
    w <- make_gradient_rast(r)
    ref <- gm_moment_of_inertia_index(w, n_bins = 100000)
    same <- gm_moment_of_inertia_index(w, n_bins = 5000)  # both above n_valid (121 cells)
    expect_equal(same$J_ref, ref$J_ref, tolerance = 1e-12)
    expect_equal(same$index, ref$index, tolerance = 1e-12)
})

test_that("forcing the binned path (n_bins < n_valid_cells) stays close to the exact reference", {
    r <- make_square_mask()
    w <- make_gradient_rast(r)
    exact <- gm_moment_of_inertia_index(w, n_bins = 100000)
    binned <- gm_moment_of_inertia_index(w, n_bins = 10)  # forces binning
    expect_true(is.finite(binned$index))
    expect_equal(binned$J_ref, exact$J_ref, tolerance = 0.05)
})

test_that("a constant weight raster does not break the binned path's degenerate range case", {
    r <- make_square_mask()
    unweighted <- gm_moment_of_inertia_index(r)
    binned_constant <- gm_moment_of_inertia_index(r, weighted = FALSE, n_bins = 5)
    expect_equal(binned_constant$index, unweighted$index, tolerance = 1e-10)
})
