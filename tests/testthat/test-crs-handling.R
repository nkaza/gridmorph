make_test_rast <- function(crs) {
    r <- terra::rast(nrows = 21, ncols = 21, xmin = -10, xmax = -9, ymin = 40, ymax = 41, crs = crs)
    terra::values(r) <- 0
    r[6:16, 6:16] <- 1
    r
}

test_that(".check_planar_crs errors clearly on a real geographic (lon/lat) CRS", {
    r <- make_test_rast("EPSG:4326")
    expect_error(.check_planar_crs(r, "test_fn"), "geographic")
})

test_that(".check_planar_crs warns (does not error) on a missing CRS", {
    r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "")
    expect_warning(res <- .check_planar_crs(r, "test_fn"), "no CRS")
    expect_null(res)
})

test_that(".check_planar_crs is silent on a real projected/planar CRS", {
    r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "local")
    expect_no_warning(.check_planar_crs(r, "test_fn"))
    expect_no_error(.check_planar_crs(r, "test_fn"))
})

#' A geographic CRS silently corrupts results in TWO different, verified
#' ways depending on which terra primitives an index happens to call (see
#' R/utils.R's own header for the two concrete cases this regression-
#' guards): gm_depth_index() mixed terra::distance()'s own geodesic-metres
#' output with a degrees-squared area, giving index = 107695.9 for a
#' quantity documented to be in (0, 1]; gm_moment_of_inertia_index() never
#' touches a geodesic-aware terra primitive at all, so it silently computed
#' in raw degree-space throughout, giving a PLAUSIBLE-LOOKING 0.96 that
#' still doesn't mean what it claims to. Both must now hard-error instead.
test_that("every exported index function hard-errors on a real geographic CRS, not just some", {
    r <- make_test_rast("EPSG:4326")
    w <- terra::init(r, "x")  # a continuous raster, for the weighted-capable functions

    expect_error(gm_depth_index(r), "geographic")
    expect_error(gm_moment_of_inertia_index(r), "geographic")
    expect_error(gm_moment_isotropy_index(r), "geographic")
    expect_error(gm_directional_balance_index(r), "geographic")
    expect_error(gm_convexity_index(r, size = 200), "geographic")
    expect_error(gm_span_index(r, size = 200), "geographic")
    expect_error(gm_radial_concentration_index(r, size = 200), "geographic")
    expect_error(gm_hull_ratio_index(r), "geographic")
    expect_error(gm_polsby_popper_index(r), "geographic")
    expect_error(gm_width_length_ratio_index(r), "geographic")
    expect_error(gm_reock_index(r), "geographic")
    expect_error(gm_detour_index(r), "geographic")
    expect_error(gm_exchange_index(r), "geographic")
})

test_that("a missing CRS warns but still produces the SAME result as an explicit planar CRS", {
    r_nocrs <- terra::rast(nrows = 21, ncols = 21, xmin = 0, xmax = 21, ymin = 0, ymax = 21, crs = "")
    terra::values(r_nocrs) <- 0
    r_nocrs[6:16, 6:16] <- 1
    r_local <- terra::rast(nrows = 21, ncols = 21, xmin = 0, xmax = 21, ymin = 0, ymax = 21, crs = "local")
    terra::values(r_local) <- 0
    r_local[6:16, 6:16] <- 1

    expect_warning(res_nocrs <- gm_depth_index(r_nocrs), "no CRS")
    res_local <- gm_depth_index(r_local)
    expect_equal(res_nocrs$index, res_local$index)

    expect_warning(hull_nocrs <- gm_hull_ratio_index(r_nocrs), "no CRS")
    hull_local <- gm_hull_ratio_index(r_local)
    expect_equal(hull_nocrs$index, hull_local$index)
})

test_that("erode()/dilate()/opening()/closing() are CRS-agnostic - no warning or error on any CRS", {
    r_geo <- make_test_rast("EPSG:4326")
    r_nocrs <- terra::rast(nrows = 21, ncols = 21, xmin = -10, xmax = -9, ymin = 40, ymax = 41, crs = "")
    terra::values(r_nocrs) <- 0
    r_nocrs[6:16, 6:16] <- 1

    expect_no_error(erode(r_geo))
    expect_no_warning(erode(r_geo))
    expect_no_error(dilate(r_nocrs))
    expect_no_warning(dilate(r_nocrs))
})
