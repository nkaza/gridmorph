make_disk_mask <- function(n = 61, radius = 20) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    cc <- terra::init(r, "x") - n / 2
    cy <- terra::init(r, "y") - n / 2
    terra::ifel(sqrt(cc^2 + cy^2) <= radius, 1, 0)
}

test_that("a disk-like raster scores close to 1", {
    res <- gm_radial_concentration_index(make_disk_mask(), size = 3000, seed = 1)
    expect_true(is.finite(res$index))
    expect_gt(res$index, 0.99)
})

test_that("an elongated rectangle scores much lower than a disk", {
    r <- terra::rast(nrows = 61, ncols = 61, xmin = 0, xmax = 61, ymin = 0, ymax = 61, crs = "local")
    terra::values(r) <- 0
    r[28:33, 5:55] <- 1
    res <- gm_radial_concentration_index(r, size = 3000, seed = 1)
    expect_lt(res$index, 0.7)
})

test_that("weighted = FALSE on a continuous raster exactly reproduces the plain binary-mask index", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[10:30, 10:30] <- 1
    w <- terra::ifel(r == 1, terra::init(r, "x") + terra::init(r, "y") + 1, NA)

    plain <- gm_radial_concentration_index(r, size = 4000, seed = 1)
    ignored <- gm_radial_concentration_index(w, weighted = FALSE, size = 4000, seed = 1)
    weighted <- gm_radial_concentration_index(w, weighted = TRUE, size = 4000, seed = 1)

    # both routes take the same "degenerate" (constant-density) sampling
    # path over the same underlying shape - bit-identical, not just
    # statistically close (see test-span-index.R's identical comment)
    expect_identical(ignored$D1_ref, plain$D1_ref)
    expect_identical(ignored$D1, plain$D1)
    expect_identical(ignored$index, plain$index)
    expect_false(isTRUE(all.equal(weighted$index, ignored$index, tolerance = 1e-6)))
})

test_that("multi-part shape's geometric median lands between the two parts, not inside either", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[5:10, 5:10] <- 1
    r[30:35, 30:35] <- 1
    res <- gm_radial_concentration_index(r, size = 3000, seed = 1)
    expect_gt(res$center["x"], 4)
    expect_lt(res$center["x"], 36)
})

test_that("size exceeding the memory-derived ceiling hard-stops with a clear error", {
    r <- make_disk_mask()
    expect_error(gm_radial_concentration_index(r, size = 1e15), "size.*exceeds")
})

test_that("no valid cells gives NA with a warning, not an error", {
    r <- make_disk_mask()
    terra::values(r) <- 0
    expect_warning(res <- gm_radial_concentration_index(r), "No valid cells")
    expect_true(is.na(res$index))
})
