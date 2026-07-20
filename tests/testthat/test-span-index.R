make_disk_mask <- function(n = 61, radius = 20) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    cc <- terra::init(r, "x") - n / 2
    cy <- terra::init(r, "y") - n / 2
    terra::ifel(sqrt(cc^2 + cy^2) <= radius, 1, 0)
}

test_that("a disk-like raster scores close to 1", {
    res <- gm_span_index(make_disk_mask(), size = 3000, seed = 1)
    expect_true(is.finite(res$index))
    expect_gt(res$index, 0.99)
})

test_that("an elongated rectangle scores much lower than a disk", {
    r <- terra::rast(nrows = 61, ncols = 61, xmin = 0, xmax = 61, ymin = 0, ymax = 61, crs = "local")
    terra::values(r) <- 0
    r[28:33, 5:55] <- 1
    res <- gm_span_index(r, size = 3000, seed = 1)
    expect_lt(res$index, 0.7)
})

test_that("weighted = FALSE on a continuous raster exactly reproduces the plain binary-mask index", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[10:30, 10:30] <- 1
    w <- terra::ifel(r == 1, terra::init(r, "x") + terra::init(r, "y") + 1, NA)

    plain <- gm_span_index(r, size = 4000, seed = 1)
    ignored <- gm_span_index(w, weighted = FALSE, size = 4000, seed = 1)
    weighted <- gm_span_index(w, weighted = TRUE, size = 4000, seed = 1)

    # both routes take the same "degenerate" (constant-density) sampling
    # path over the same underlying shape, so this is bit-identical, not
    # just statistically close - see R/utils.R's .sample_valid_points()
    expect_identical(ignored$D_ref, plain$D_ref)
    expect_identical(ignored$D, plain$D)
    expect_identical(ignored$index, plain$index)
    # and the ACTUAL weighted computation differs (the gradient is genuinely non-uniform)
    expect_false(isTRUE(all.equal(weighted$index, ignored$index, tolerance = 1e-6)))
})

test_that("D_ref depends only on the weight histogram, not where it's placed", {
    r <- terra::rast(nrows = 61, ncols = 61, xmin = 0, xmax = 61, ymin = 0, ymax = 61, crs = "local")
    terra::values(r) <- 0
    r[8:53, 8:53] <- 1
    valid <- !is.na(r) & (r == 1)
    valid_v <- as.vector(terra::values(valid))

    cc <- terra::init(r, "x") - 30.5
    cy <- terra::init(r, "y") - 30.5
    dist_c <- sqrt(cc^2 + cy^2)
    d_valid <- as.vector(terra::values(dist_c))[valid_v]
    ranks <- rank(d_valid, ties.method = "first")
    n <- length(ranks)

    w_center_full <- rep(NA_real_, terra::ncell(r))
    w_edge_full <- rep(NA_real_, terra::ncell(r))
    w_center_full[valid_v] <- n - ranks + 1  # high weight near centre
    w_edge_full[valid_v] <- ranks            # high weight near edge

    w_center <- r; terra::values(w_center) <- w_center_full
    w_edge <- r; terra::values(w_edge) <- w_edge_full

    res_center <- gm_span_index(w_center, size = 4000, seed = 1)
    res_edge <- gm_span_index(w_edge, size = 4000, seed = 1)

    expect_equal(res_center$D_ref, res_edge$D_ref, tolerance = 1e-8)
    # concentrating mass near the centre should score closer to 1 than
    # concentrating it near the edge (Riesz rearrangement direction)
    expect_gt(res_center$index, res_edge$index)
})

test_that("size exceeding the memory-derived ceiling hard-stops with a clear error", {
    r <- make_disk_mask()
    expect_error(gm_span_index(r, size = 1e15), "size.*exceeds")
})

test_that("no valid cells gives NA with a warning, not an error", {
    r <- make_disk_mask()
    terra::values(r) <- 0
    expect_warning(res <- gm_span_index(r), "No valid cells")
    expect_true(is.na(res$index))
})
