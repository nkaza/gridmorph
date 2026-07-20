make_square_mask <- function(n = 41, inset = 10) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    terra::values(r) <- 0
    r[(inset + 1):(n - inset), (inset + 1):(n - inset)] <- 1
    r
}

test_that("a solid square scores (close to) 1, fully convex", {
    r <- make_square_mask()
    res <- gm_convexity_index(r, size = 2000, seed = 1)
    expect_true(is.finite(res$index))
    expect_gt(res$index, 0.99)
    expect_lte(res$index, 1)
})

test_that("two disjoint blobs score much lower than one solid square", {
    r <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r) <- 0
    r[5:10, 5:10] <- 1
    r[30:35, 30:35] <- 1
    res <- gm_convexity_index(r, size = 2000, seed = 1)
    sq <- gm_convexity_index(make_square_mask(), size = 2000, seed = 1)
    expect_lt(res$index, sq$index - 0.2)
})

test_that("a hole through the middle lowers convexity relative to no hole", {
    r <- make_square_mask()
    r_hole <- r
    r_hole[16:26, 16:26] <- NA  # hole
    with_hole <- gm_convexity_index(r_hole, size = 3000, seed = 1)
    no_hole <- gm_convexity_index(r, size = 3000, seed = 1)
    expect_lt(with_hole$index, no_hole$index)
})

test_that("reproducible with the same seed, different across seeds", {
    r <- make_square_mask(n = 61, inset = 5)  # a bit irregular-shaped to give seed-dependent variance
    r[5:15, 45:55] <- 0  # notch a corner to introduce concavity
    a <- gm_convexity_index(r, size = 500, seed = 42)$index
    b <- gm_convexity_index(r, size = 500, seed = 42)$index
    c <- gm_convexity_index(r, size = 500, seed = 43)$index
    expect_identical(a, b)
    expect_false(isTRUE(all.equal(a, c)))
})

test_that("size exceeding the memory-derived ceiling hard-stops with a clear error", {
    r <- make_square_mask()
    expect_error(gm_convexity_index(r, size = 1e15), "size.*exceeds")
})

test_that("a single valid cell is vacuously convex (index = 1) without sampling", {
    r <- make_square_mask()
    terra::values(r) <- 0
    r[20, 20] <- 1
    res <- gm_convexity_index(r)
    expect_equal(res$index, 1)
    expect_equal(res$n_valid_cells, 1L)
})

test_that("no valid cells gives NA with a warning, not an error", {
    r <- make_square_mask()
    terra::values(r) <- 0
    expect_warning(res <- gm_convexity_index(r), "No valid cells")
    expect_true(is.na(res$index))
})
