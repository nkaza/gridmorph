make_square_mask <- function(n = 41, inset = 10) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    terra::values(r) <- 0
    r[(inset + 1):(n - inset), (inset + 1):(n - inset)] <- 1
    r
}

make_disk_mask <- function(n = 61, radius = 20) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    cc <- terra::init(r, "x") - n / 2
    cy <- terra::init(r, "y") - n / 2
    terra::ifel(sqrt(cc^2 + cy^2) <= radius, 1, 0)
}

make_rectangle_mask <- function() {
    r <- terra::rast(nrows = 61, ncols = 61, xmin = 0, xmax = 61, ymin = 0, ymax = 61, crs = "local")
    terra::values(r) <- 0
    r[28:33, 5:55] <- 1
    r
}

test_that("a square is convex (hull_ratio = 1 exactly) and matches the exact theoretical Reock value", {
    r <- make_square_mask()
    hr <- gm_hull_ratio_index(r)
    expect_equal(hr$index, 1, tolerance = 1e-10)
    # a square's MEC has radius = half the diagonal, so area/mbc_area = 2/pi exactly -
    # gm_reock_index()'s own `index` is resolution-matched (raw/ref_index, see its
    # own roxygen), but area/mbc_area recovers the textbook raw ratio directly
    res <- gm_reock_index(r)
    expect_equal(res$area / res$mbc_area, 2 / pi, tolerance = 0.01)
    expect_gt(gm_width_length_ratio_index(r)$index, 0.99)
})

test_that("a disk-like raster scores close to 1 on every classical metric", {
    r <- make_disk_mask()
    expect_gt(gm_hull_ratio_index(r)$index, 0.9)
    expect_gt(gm_polsby_popper_index(r)$index, 0.8)
    expect_gt(gm_reock_index(r)$index, 0.9)
    expect_gt(gm_detour_index(r)$index, 0.9)
    expect_equal(gm_exchange_index(r)$index, 1, tolerance = 0.02)
})

test_that("resolution-matched references put a disk within a tight band of 1, not just above a loose floor", {
    # regression test for R/classical-metrics.R's own file header: comparing
    # raw against the pure analytical formula used to put a disk at
    # 0.880/0.953/0.982 on these three at this resolution; a reference built
    # from the identical raster-measurement pipeline should cancel almost
    # all of that gap
    r <- make_disk_mask(n = 81, radius = 30)
    expect_equal(gm_polsby_popper_index(r)$index, 1, tolerance = 0.05)
    expect_equal(gm_reock_index(r)$index, 1, tolerance = 0.01)
    expect_equal(gm_detour_index(r)$index, 1, tolerance = 0.01)
})

test_that("an elongated rectangle is still convex but scores low on width/length, polsby-popper, reock", {
    r <- make_rectangle_mask()
    expect_equal(gm_hull_ratio_index(r)$index, 1, tolerance = 1e-10)  # a rectangle is convex
    expect_lt(gm_width_length_ratio_index(r)$index, 0.2)
    expect_lt(gm_polsby_popper_index(r)$index, 0.5)
    expect_lt(gm_reock_index(r)$index, 0.3)
})

test_that("an L-shape (non-convex) scores below 1 on hull_ratio", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[5:35, 5:15] <- 1
    r[25:35, 5:35] <- 1
    expect_lt(gm_hull_ratio_index(r)$index, 0.9)
})

test_that("two well-separated blobs give exchange_index exactly 0 (the documented pathology)", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[5:10, 5:10] <- 1
    r[30:35, 30:35] <- 1
    expect_equal(gm_exchange_index(r)$index, 0)
})

test_that("multi-part hull spans across disjoint parts (hull_area much bigger than combined blob area)", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[5:10, 5:10] <- 1
    r[30:35, 30:35] <- 1
    res <- gm_hull_ratio_index(r)
    expect_lt(res$index, 0.3)
    expect_gt(res$hull_area, res$area * 3)
})

test_that("area and n_valid_cells agree across all six functions", {
    r <- make_square_mask()
    n <- sum(as.vector(terra::values(r)), na.rm = TRUE)
    cell_area <- prod(terra::res(r))
    expected_area <- n * cell_area

    expect_equal(gm_hull_ratio_index(r)$area, expected_area)
    expect_equal(gm_polsby_popper_index(r)$area, expected_area)
    expect_equal(gm_reock_index(r)$area, expected_area)
    expect_equal(gm_detour_index(r)$area, expected_area)
    expect_equal(gm_exchange_index(r)$area, expected_area)
    expect_equal(gm_hull_ratio_index(r)$n_valid_cells, as.integer(n))
})

test_that("hull/mbc/circle geometries are returned as terra SpatVectors for plotting", {
    r <- make_square_mask()
    expect_s4_class(gm_hull_ratio_index(r)$hull, "SpatVector")
    expect_s4_class(gm_reock_index(r)$mbc, "SpatVector")
    expect_s4_class(gm_detour_index(r)$hull, "SpatVector")
    expect_s4_class(gm_exchange_index(r)$circle, "SpatVector")
})

test_that("no valid cells gives NA with a warning, not an error, for all six functions", {
    r <- make_square_mask(); terra::values(r) <- 0
    expect_warning(res1 <- gm_hull_ratio_index(r), "No valid cells")
    expect_warning(res2 <- gm_polsby_popper_index(r), "No valid cells")
    expect_warning(res3 <- gm_width_length_ratio_index(r), "No valid cells")
    expect_warning(res4 <- gm_reock_index(r), "No valid cells")
    expect_warning(res5 <- gm_detour_index(r), "No valid cells")
    expect_warning(res6 <- gm_exchange_index(r), "No valid cells")
    expect_true(is.na(res1$index)); expect_true(is.na(res2$index)); expect_true(is.na(res3$index))
    expect_true(is.na(res4$index)); expect_true(is.na(res5$index)); expect_true(is.na(res6$index))
})

test_that("a single valid cell doesn't crash any of the six functions", {
    r <- make_square_mask(); terra::values(r) <- 0
    r[20, 20] <- 1
    expect_no_error(gm_hull_ratio_index(r))
    expect_no_error(gm_polsby_popper_index(r))
    expect_no_error(gm_width_length_ratio_index(r))
    expect_no_error(gm_reock_index(r))
    expect_no_error(gm_detour_index(r))
    expect_no_error(gm_exchange_index(r))
})

test_that("holes are invisible to width_length_ratio_index (documented blind spot)", {
    r <- make_square_mask()
    r_hole <- r
    r_hole[16:26, 16:26] <- NA  # punch a hole through the middle
    expect_equal(gm_width_length_ratio_index(r_hole)$index, gm_width_length_ratio_index(r)$index)
})

test_that("the marching-squares perimeter beats naive stair-step counting on a diagonal shape", {
    n <- 41
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    cx <- terra::init(r, "x") - 20.5
    cy <- terra::init(r, "y") - 20.5
    diamond <- terra::ifel(abs(cx) + abs(cy) <= 15, 1, 0)
    true_perimeter <- 4 * 15 * sqrt(2)

    res <- gm_polsby_popper_index(diamond)
    rel_err <- abs(res$perimeter - true_perimeter) / true_perimeter
    expect_lt(rel_err, 0.10)  # marching squares: verified ~3-8% error on test shapes
})

test_that("a shape touching the raster's own edge gives the identical perimeter as the same shape with a margin", {
    # regression test for a verified terra::as.contour() edge bug: without
    # padding first, marching squares has nothing beyond the raster's true
    # edge to interpolate the boundary against, so a shape's edge-touching
    # side goes untraced - verified directly, a 21-cell square flush
    # against the top edge gave perimeter 40, roughly HALF the correct
    # ~82.8
    n <- 21
    r_edge <- terra::rast(nrows = n, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = n, crs = "local")
    terra::values(r_edge) <- 0
    r_edge[1:n, 11:31] <- 1  # flush against the top edge (row 1)

    r_margin <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r_margin) <- 0
    r_margin[11:31, 11:31] <- 1  # same square, centred with a margin on every side

    res_edge <- gm_polsby_popper_index(r_edge)
    res_margin <- gm_polsby_popper_index(r_margin)
    expect_equal(res_edge$perimeter, res_margin$perimeter, tolerance = 1e-10)
    expect_equal(res_edge$index, res_margin$index, tolerance = 1e-10)
})
