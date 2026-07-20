make_disk_mask <- function(n = 61, radius = 20) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    cc <- terra::init(r, "x") - n / 2
    cy <- terra::init(r, "y") - n / 2
    terra::ifel(sqrt(cc^2 + cy^2) <= radius, 1, 0)
}

## -- gm_geodesic_span_index() ---------------------------------------------

test_that("a disk-like raster scores close to 1 on geodesic_span", {
    res <- gm_geodesic_span_index(make_disk_mask(), n_points = 150, seed = 1)
    expect_true(is.finite(res$index))
    expect_equal(res$index, 1, tolerance = 0.1)
})

test_that("an elongated rectangle scores much lower than a disk on geodesic_span", {
    r <- terra::rast(nrows = 61, ncols = 61, xmin = 0, xmax = 61, ymin = 0, ymax = 61, crs = "local")
    terra::values(r) <- 0
    r[28:33, 5:55] <- 1
    res <- gm_geodesic_span_index(r, n_points = 100, seed = 1)
    expect_lt(res$index, 0.75)
})

test_that("weighted = FALSE on a continuous raster exactly reproduces the plain binary-mask geodesic_span index", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[10:30, 10:30] <- 1
    w <- terra::ifel(r == 1, terra::init(r, "x") + terra::init(r, "y") + 1, NA)

    plain <- gm_geodesic_span_index(r, n_points = 30, seed = 1)
    ignored <- gm_geodesic_span_index(w, weighted = FALSE, n_points = 30, seed = 1)
    weighted <- gm_geodesic_span_index(w, weighted = TRUE, n_points = 30, seed = 1)

    # both routes draw from the same degenerate (constant-density) code
    # path over the same shape - bit-identical, not just statistically
    # close, the same way gm_span_index()'s own equivalent test is
    expect_identical(ignored$D, plain$D)
    expect_identical(ignored$D_ref, plain$D_ref)
    expect_identical(ignored$index, plain$index)
    # and the ACTUAL weighted computation differs (the gradient is genuinely non-uniform)
    expect_false(isTRUE(all.equal(weighted$index, ignored$index, tolerance = 1e-6)))
})

test_that("a shape touching the raster's own edge gives the identical geodesic_span D as the same shape with a margin", {
    n <- 21
    r_edge <- terra::rast(nrows = n, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = n, crs = "local")
    terra::values(r_edge) <- 0
    r_edge[1:n, 11:31] <- 1

    r_margin <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r_margin) <- 0
    r_margin[11:31, 11:31] <- 1

    res_edge <- gm_geodesic_span_index(r_edge, n_points = 20, seed = 1)
    res_margin <- gm_geodesic_span_index(r_margin, n_points = 20, seed = 1)
    expect_equal(res_edge$D, res_margin$D, tolerance = 1e-10)
    expect_equal(res_edge$index, res_margin$index, tolerance = 1e-10)
})

test_that("n_points exceeding the memory/time-derived ceiling hard-stops with a clear error, geodesic_span", {
    r <- make_disk_mask()
    expect_error(gm_geodesic_span_index(r, n_points = 1e12), "n_points.*exceeds")
})

test_that("no valid cells gives NA with a warning, not an error, geodesic_span", {
    r <- make_disk_mask()
    terra::values(r) <- 0
    expect_warning(res <- gm_geodesic_span_index(r), "No valid cells")
    expect_true(is.na(res$index))
})

test_that("a single valid cell gives NA with a warning, not an error, geodesic_span", {
    r <- make_disk_mask()
    terra::values(r) <- 0
    r[30, 30] <- 1
    expect_warning(res <- gm_geodesic_span_index(r), "one valid cell")
    expect_true(is.na(res$index))
    expect_equal(res$n_valid_cells, 1L)
})

test_that("a hole punched into a disk does not crash geodesic_span and lowers the index", {
    r <- make_disk_mask(n = 61, radius = 25)
    r_hole <- r
    r_hole[27:33, 27:33] <- NA
    plain <- gm_geodesic_span_index(r, n_points = 80, seed = 1)
    holed <- gm_geodesic_span_index(r_hole, n_points = 80, seed = 1)
    expect_true(is.finite(holed$index))
    expect_lt(holed$n_valid_cells, plain$n_valid_cells)
})

test_that("two disjoint parts warn with a specific reason and return NA, not a bare NaN, for geodesic_span", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- NA
    r[5:10, 5:10] <- 1
    r[30:35, 30:35] <- 1
    # caught by the cheap up-front .is_connected() check, before any
    # gridDist() sampling runs at all - so D_ref is NA too (never built),
    # not finite the way it would be if only a sampled PAIR had failed
    expect_warning(res <- gm_geodesic_span_index(r, n_points = 10, seed = 1), "more than one connected part")
    expect_true(is.na(res$index))
    expect_true(is.na(res$D))
    expect_true(is.na(res$D_ref))
    expect_false(is.nan(res$index))
    expect_false(is.nan(res$D))
    expect_false(is.nan(res$D_ref))
})

## -- gm_geodesic_chord_index() ---------------------------------------------

test_that("a disk-like raster scores close to 1 on geodesic_chord", {
    res <- gm_geodesic_chord_index(make_disk_mask(), n_points = 100, seed = 1)
    expect_true(is.finite(res$index))
    expect_equal(res$index, 1, tolerance = 0.08)
})

test_that("gm_geodesic_chord_index() has no weighted argument at all", {
    expect_false("weighted" %in% names(formals(gm_geodesic_chord_index)))
    r <- make_disk_mask()
    expect_error(gm_geodesic_chord_index(r, weighted = TRUE), "unused argument")
})

test_that("a shape touching the raster's own edge gives the identical geodesic_chord D as the same shape with a margin", {
    n <- 21
    r_edge <- terra::rast(nrows = n, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = n, crs = "local")
    terra::values(r_edge) <- 0
    r_edge[1:n, 11:31] <- 1

    r_margin <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r_margin) <- 0
    r_margin[11:31, 11:31] <- 1

    res_edge <- gm_geodesic_chord_index(r_edge, n_points = 20, seed = 1)
    res_margin <- gm_geodesic_chord_index(r_margin, n_points = 20, seed = 1)
    expect_equal(res_edge$D, res_margin$D, tolerance = 1e-10)
    expect_equal(res_edge$index, res_margin$index, tolerance = 1e-10)
})

test_that("n_points exceeding the number of boundary cells warns and falls back to sampling all of them, geodesic_chord", {
    r <- make_disk_mask(n = 21, radius = 5)
    n_bnd <- sum(as.vector(terra::values(gridmorph:::.boundary_cells(gridmorph:::.valid_cells(r)))))
    expect_warning(res <- gm_geodesic_chord_index(r, n_points = n_bnd + 50, seed = 1), "exceeds the number of boundary cells")
    expect_equal(res$n_boundary_cells, n_bnd)
})

test_that("n_points exceeding the memory/time-derived ceiling hard-stops with a clear error, geodesic_chord", {
    r <- make_disk_mask()
    expect_error(gm_geodesic_chord_index(r, n_points = 1e12), "n_points.*exceeds")
})

test_that("no valid cells gives NA with a warning, not an error, geodesic_chord", {
    r <- make_disk_mask()
    terra::values(r) <- 0
    expect_warning(res <- gm_geodesic_chord_index(r), "No valid cells")
    expect_true(is.na(res$index))
})

test_that("an elongated rectangle scores lower than a disk on geodesic_chord", {
    r <- terra::rast(nrows = 61, ncols = 61, xmin = 0, xmax = 61, ymin = 0, ymax = 61, crs = "local")
    terra::values(r) <- 0
    r[28:33, 5:55] <- 1
    res <- gm_geodesic_chord_index(r, n_points = 80, seed = 1)
    expect_lt(res$index, 0.9)
})

test_that("two disjoint parts warn with a specific reason and return NA, not a bare NaN, for geodesic_chord", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- NA
    r[5:10, 5:10] <- 1
    r[30:35, 30:35] <- 1
    # caught by the cheap up-front .is_connected() check - D_ref is NA too
    # (never built), n_boundary_cells is still real (computed before the
    # connectivity check, cheap either way)
    expect_warning(res <- gm_geodesic_chord_index(r, n_points = 10, seed = 1), "more than one connected part")
    expect_true(is.na(res$index))
    expect_true(is.na(res$D))
    expect_true(is.na(res$D_ref))
    expect_false(is.nan(res$index))
    expect_false(is.nan(res$D))
    expect_false(is.nan(res$D_ref))
    expect_gt(res$n_boundary_cells, 0)
})

test_that(".is_connected() matches gridDist()'s own 8-connectivity, not 4-connectivity", {
    # two 3x3 blocks touching ONLY at a single shared corner (the diagonal
    # pinch between cell (5,5) and cell (6,6)) - reachable via gridDist()'s
    # own diagonal moves (verified in R/geodesic-index.R's own file
    # header), so this must NOT warn or return NA. Multi-cell blocks
    # (rather than two lone cells) avoid an unrelated small-sample
    # duplicate-point edge case in .sample_valid_points() at n_valid = 2.
    r <- terra::rast(nrows = 15, ncols = 15, xmin = 0, xmax = 15, ymin = 0, ymax = 15, crs = "local")
    terra::values(r) <- NA
    r[3:5, 3:5] <- 1
    r[6:8, 6:8] <- 1
    expect_no_warning(res <- gm_geodesic_span_index(r, n_points = 8, seed = 1))
    expect_true(is.finite(res$index))
})

test_that("a connectivity-disqualified shape skips the expensive sampling loop entirely", {
    # regression guard for the up-front .is_connected() check actually
    # short-circuiting - if it didn't, this would still run K sequential
    # whole-raster gridDist() calls before discovering the same answer
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- NA
    r[5:10, 5:10] <- 1
    r[30:35, 30:35] <- 1
    called <- 0L
    testthat::local_mocked_bindings(
        .geodesic_pairwise_mean = function(...) { called <<- called + 1L; list(D = NA_real_, K = 0L, n_unreachable_pairs = 0L) },
        .package = "gridmorph"
    )
    suppressWarnings(gm_geodesic_span_index(r, n_points = 10, seed = 1))
    expect_equal(called, 0L)
})

## -- shared helper -----------------------------------------------------

test_that(".geodesic_pairwise_mean() de-duplicates repeated cell numbers before pairing", {
    r <- make_disk_mask(n = 21, radius = 8)
    valid <- gridmorph:::.valid_cells(r)
    cells <- which(as.vector(terra::values(valid)))[1:5]
    with_dupe <- gridmorph:::.geodesic_pairwise_mean(valid, c(cells, cells[1]))
    without_dupe <- gridmorph:::.geodesic_pairwise_mean(valid, cells)
    expect_equal(with_dupe$K, without_dupe$K)
    expect_equal(with_dupe$D, without_dupe$D)
})

test_that(".geodesic_pairwise_mean() returns NA with K < 2 for fewer than two unique cells", {
    r <- make_disk_mask(n = 21, radius = 8)
    valid <- gridmorph:::.valid_cells(r)
    cell <- which(as.vector(terra::values(valid)))[1]
    res <- gridmorph:::.geodesic_pairwise_mean(valid, c(cell, cell))
    expect_true(is.na(res$D))
    expect_equal(res$K, 1L)
})
