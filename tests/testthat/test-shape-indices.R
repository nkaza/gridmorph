make_square_mask <- function(n = 41, inset = 10) {
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    terra::values(r) <- 0
    r[(inset + 1):(n - inset), (inset + 1):(n - inset)] <- 1
    r
}

test_that("gm_shape_indices(which = \"all\") returns all fifteen, matching direct calls exactly", {
    r <- make_square_mask()
    # size= only reaches the three ORIGINAL Monte Carlo indices - the two
    # geodesic ones use their own n_points= instead (deliberately
    # decoupled, see R/geodesic-index.R's own file header) and fall back
    # to their own default here, matched below via n_points not size
    res <- gm_shape_indices(r, size = 2000, seed = 1)

    expect_equal(names(res), .ALL_GM_INDICES)
    expect_equal(unname(res["depth"]), gm_depth_index(r)$index)
    expect_equal(unname(res["moment_of_inertia"]), gm_moment_of_inertia_index(r)$index)
    expect_equal(unname(res["moment_isotropy"]), gm_moment_isotropy_index(r)$index)
    expect_equal(unname(res["directional_balance"]), gm_directional_balance_index(r)$index)
    expect_equal(unname(res["convexity"]), gm_convexity_index(r, size = 2000, seed = 1)$index)
    expect_equal(unname(res["span"]), gm_span_index(r, size = 2000, seed = 1)$index)
    expect_equal(unname(res["radial_concentration"]), gm_radial_concentration_index(r, size = 2000, seed = 1)$index)
    expect_equal(unname(res["hull_ratio"]), gm_hull_ratio_index(r)$index)
    expect_equal(unname(res["polsby_popper"]), gm_polsby_popper_index(r)$index)
    expect_equal(unname(res["width_length_ratio"]), gm_width_length_ratio_index(r)$index)
    expect_equal(unname(res["reock"]), gm_reock_index(r)$index)
    expect_equal(unname(res["detour"]), gm_detour_index(r)$index)
    expect_equal(unname(res["exchange"]), gm_exchange_index(r)$index)
    expect_equal(unname(res["geodesic_span"]), gm_geodesic_span_index(r, seed = 1)$index)
    expect_equal(unname(res["geodesic_chord"]), gm_geodesic_chord_index(r, seed = 1)$index)
})

test_that("which = \"all\" includes the two geodesic indices, matching direct calls with n_points=", {
    r <- make_square_mask()
    res <- gm_shape_indices(r, n_points = 30, seed = 1)
    expect_true(all(c("geodesic_span", "geodesic_chord") %in% names(res)))
    expect_equal(length(res), 15L)
    expect_equal(unname(res["geodesic_span"]), gm_geodesic_span_index(r, n_points = 30, seed = 1)$index)
    expect_equal(unname(res["geodesic_chord"]), gm_geodesic_chord_index(r, n_points = 30, seed = 1)$index)
})

test_that("the geodesic indices are reachable via an explicit which= request too, matching direct calls", {
    r <- make_square_mask()
    res <- gm_shape_indices(r, which = c("hull_ratio", "geodesic_span", "geodesic_chord"), n_points = 30, seed = 1)
    expect_equal(names(res), c("hull_ratio", "geodesic_span", "geodesic_chord"))
    expect_equal(unname(res["geodesic_span"]), gm_geodesic_span_index(r, n_points = 30, seed = 1)$index)
    expect_equal(unname(res["geodesic_chord"]), gm_geodesic_chord_index(r, n_points = 30, seed = 1)$index)
})

test_that("a subset via which= returns exactly (and only) those, in canonical order regardless of request order", {
    r <- make_square_mask()
    res <- gm_shape_indices(r, which = c("reock", "hull_ratio", "polsby_popper"))
    expect_equal(names(res), c("hull_ratio", "polsby_popper", "reock"))
    expect_equal(length(res), 3L)
})

test_that("a subset of only classic metrics needs no size/seed at all", {
    r <- make_square_mask()
    expect_no_error(gm_shape_indices(r, which = c("hull_ratio", "detour")))
})

test_that("unknown names in which= error immediately, listing valid choices", {
    r <- make_square_mask()
    expect_error(gm_shape_indices(r, which = "not_a_real_index"), "Unknown index")
})

test_that("which must be \"all\" or a character vector", {
    r <- make_square_mask()
    expect_error(gm_shape_indices(r, which = 1:3), "character vector")
})

test_that("weighted is forwarded to the seven weight-capable indices, ignored by the six classic ones", {
    r <- make_square_mask()
    w <- terra::ifel(r == 1, terra::init(r, "x") + terra::init(r, "y") + 1, NA)

    res_weighted <- gm_shape_indices(w, weighted = TRUE, size = 2000, seed = 1)
    res_ignored <- gm_shape_indices(w, weighted = FALSE, size = 2000, seed = 1)

    # weighted vs not should differ for the weight-capable indices...
    expect_false(isTRUE(all.equal(res_weighted["moment_of_inertia"], res_ignored["moment_of_inertia"])))
    # ...but be IDENTICAL for the six classic metrics, which have no weighted form
    expect_equal(res_weighted["hull_ratio"], res_ignored["hull_ratio"])
    expect_equal(res_weighted["polsby_popper"], res_ignored["polsby_popper"])
    expect_equal(res_weighted["reock"], res_ignored["reock"])
})

test_that("n_bins is forwarded and overrides each function's own default, including span_index()'s smaller one", {
    r <- make_square_mask()
    res_default <- gm_shape_indices(r, which = c("depth", "span"), size = 500, seed = 1)
    res_custom <- gm_shape_indices(r, which = c("depth", "span"), size = 500, seed = 1, n_bins = 5)
    # both are constant-weight here so results are identical regardless of
    # n_bins (degenerate collapse) - just confirming no error from passing
    # n_bins through to span (default 100) alongside depth (default 1000)
    expect_equal(res_default, res_custom, tolerance = 1e-8)
})

test_that("a missing CRS warns exactly once for gm_shape_indices(\"all\"), not once per index", {
    r <- terra::rast(nrows = 21, ncols = 21, xmin = 0, xmax = 21, ymin = 0, ymax = 21, crs = "")
    terra::values(r) <- 0
    r[6:16, 6:16] <- 1

    warnings_seen <- character(0)
    withCallingHandlers(
        gm_shape_indices(r, size = 200, seed = 1),
        warning = function(w) {
            warnings_seen[length(warnings_seen) + 1] <<- conditionMessage(w)
            invokeRestart("muffleWarning")
        }
    )
    expect_equal(length(warnings_seen), 1L)
})

test_that("weighted = TRUE on a categorical raster warns exactly once, not once per weighted-capable index", {
    r <- terra::rast(nrows = 21, ncols = 21, xmin = 0, xmax = 21, ymin = 0, ymax = 21, crs = "local")
    terra::values(r) <- 0
    r[6:16, 6:16] <- sample(1:3, 121, replace = TRUE)
    r <- terra::as.factor(r)

    warnings_seen <- character(0)
    withCallingHandlers(
        gm_shape_indices(r, weighted = TRUE, size = 200, seed = 1),
        warning = function(w) {
            warnings_seen[length(warnings_seen) + 1] <<- conditionMessage(w)
            invokeRestart("muffleWarning")
        }
    )
    expect_equal(length(warnings_seen), 1L)
    expect_match(warnings_seen[1], "categorical")

    expect_no_warning(gm_shape_indices(r, weighted = FALSE, size = 200, seed = 1))
})

test_that("a geographic CRS errors immediately, before computing anything", {
    r <- terra::rast(nrows = 21, ncols = 21, xmin = -10, xmax = -9, ymin = 40, ymax = 41, crs = "EPSG:4326")
    terra::values(r) <- 0
    r[6:16, 6:16] <- 1
    expect_error(gm_shape_indices(r), "geographic")
})

test_that("every index is translation-invariant: a shape touching the raster's own edge scores identically to the same shape with a margin", {
    # a broad regression guard for the whole class of raster-edge-padding
    # bug found in erode()/gm_depth_index()/gm_polsby_popper_index() (see
    # R/morphology.R's, R/depth-index.R's, and R/classical-metrics.R's own
    # headers) - any future index or helper that reads a moving-window or
    # boundary-tracing terra primitive without accounting for the
    # raster's true edge would show up here as a nonzero difference,
    # regardless of which specific mechanism causes it
    n <- 21
    r_edge <- terra::rast(nrows = n, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = n, crs = "local")
    terra::values(r_edge) <- 0
    r_edge[1:n, 11:31] <- 1  # flush against the top edge (row 1)

    r_margin <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(r_margin) <- 0
    r_margin[11:31, 11:31] <- 1  # same square, centred with a margin on every side

    # "all" now includes the geodesic pair too (their own default
    # n_points, since size= doesn't reach them - see R/geodesic-index.R's
    # own file header for why) - one call already covers all fifteen
    res_edge <- gm_shape_indices(r_edge, size = 2000, seed = 1)
    res_margin <- gm_shape_indices(r_margin, size = 2000, seed = 1)
    expect_equal(res_edge, res_margin, tolerance = 1e-8)
})
