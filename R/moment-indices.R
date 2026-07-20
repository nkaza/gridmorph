## =========================================================================
## Moment of inertia, moment isotropy, and directional balance indices,
## ported from shapeindices::moment_of_inertia_index()/moment_isotropy_index()/
## directional_balance_index(). All three share one mass-centroid + inertia-
## tensor core, computed via terra coordinate rasters (terra::init(x, "x"/"y"))
## and terra::global() reductions rather than pulling terra::values() into
## memory - chunk-safe for free the same way terra::distance()/patches()
## already are, per the plan's memory design.
##
## UNITS: each valid cell contributes mass = density * cell_area (cell_area
## is a constant for a raster, unlike triangle area in the vector package) -
## Ixx/Iyy/Ixy each need an explicit `* cell_area` factor the centroid
## calculation itself doesn't (it cancels in a ratio of two cell_area-scaled
## sums). Leading-order approximation: each cell's own second moment about
## its own centre (cell_area * (cell_width^2 + cell_height^2)/12) is not
## added - negligible at realistic resolutions, the raster analogue of a
## small quadrature bias already documented elsewhere in this package,
## not attempted to correct here.
##
## gm_moment_of_inertia_index(): index = J_ref/J; J_ref via the same
## concentric-rings construction as gm_depth_index()'s own weighted
## reference - identical "sort by density descending, ring boundaries from
## cumulative area" pattern. The ring construction ITSELF is shared
## (R/utils.R's `.adaptive_density_bins()` - exact per-cell sort when
## cheap, a memory-bounded K-bin histogram approximation above that
## threshold); only the final rho*(S^2-S_prev^2)/(2pi) formula stays
## separate from depth's rho*mean_r construction, since the two indices do
## different things with the same rings. ALWAYS the "weighted" formula,
## `weighted = FALSE` included (see `.mass_raster()`, R/utils.R): a
## constant-1 density collapses this ring construction to `area^2/(2*pi)`
## EXACTLY (single degenerate ring: S = area, S_prev = 0, so
## `J_ref = 1*(area^2 - 0)/(2*pi)`) - verified algebraically, not just
## assumed, so there's no separate closed-form path to maintain.
##
## gm_moment_isotropy_index(): lambda_min/lambda_max of [[Ixx,Ixy],[Ixy,Iyy]] -
## no reference shape at all, an elementary PSD-matrix fact (M is an
## integral of rank-1 PSD outer products), same as the vector package. A
## ratio of eigenvalues of the SAME matrix, so scale-invariant to whatever
## constant rho gets normalised to - `weighted = FALSE`'s constant density
## needs no special-casing here either.
##
## gm_directional_balance_index(): index = 1 - |mean(rho_norm * exp(i*theta))|,
## theta = bearing from the mass centroid. Already a normalised MEAN
## (divided by rho's own total), so also scale-invariant - same reasoning.

#' @noRd
.coord_rasters <- function(valid) {
    list(x = terra::init(valid, "x"), y = terra::init(valid, "y"))
}

#' The ONE density used everywhere in this file: `rho = (w / sum(w)) /
#' cell_area` - normalised mass divided by area, the direct raster
#' analogue of the vector package's own `rho <- w/tri_area` with
#' `w <- .normalize_weight(weight)`. `w` is always a concrete raster here
#' (either `rast` itself or a constant-1 stand-in - see `.mass_raster()`),
#' never `NULL`, so there is only ever this one code path.
#'
#' A REAL BUG lived here before this was unified, worth recording: Ixx/Iyy/Ixy
#' used to use the RAW weight raster value directly (not divided by
#' cell_area, not normalised) while J_ref's own ring construction
#' normalised the weight first - two different, inconsistent density
#' scales feeding one ratio (`J_ref/J`), caught only because a "uniform
#' weight reproduces the unweighted index exactly" test failed (0.01
#' instead of 0.96). Both Ixx/Iyy/Ixy and J_ref now route through this
#' single function so that can't happen again by construction, not just
#' by re-checking by hand.
#' @noRd
.density_raster <- function(valid, w, cell_area) {
    w <- terra::ifel(valid, w, NA)
    W <- as.numeric(terra::global(w, "sum", na.rm = TRUE)[1, 1])
    (w / W) / cell_area
}

#' @noRd
.mass_centroid_raster <- function(valid, rho, cc) {
    Wr <- as.numeric(terra::global(rho, "sum", na.rm = TRUE)[1, 1])
    Gx <- as.numeric(terra::global(cc$x * rho, "sum", na.rm = TRUE)[1, 1]) / Wr
    Gy <- as.numeric(terra::global(cc$y * rho, "sum", na.rm = TRUE)[1, 1]) / Wr
    c(x = Gx, y = Gy)
}

#' Shared core: density, mass centroid, and Ixx/Iyy/Ixy about it - `rho`
#' (not raw `w`) is what every downstream computation in this file uses,
#' including `.J_ref_raster()` below, so the same density scale appears
#' consistently in every numerator and denominator.
#' @noRd
.moment_tensor_raster <- function(valid, w, cell_area) {
    cc <- .coord_rasters(valid)
    rho <- .density_raster(valid, w, cell_area)
    G <- .mass_centroid_raster(valid, rho, cc)
    xr <- terra::ifel(valid, cc$x - G["x"], NA)
    yr <- terra::ifel(valid, cc$y - G["y"], NA)

    Ixx <- cell_area * as.numeric(terra::global(rho * yr^2, "sum", na.rm = TRUE)[1, 1])
    Iyy <- cell_area * as.numeric(terra::global(rho * xr^2, "sum", na.rm = TRUE)[1, 1])
    Ixy <- cell_area * as.numeric(terra::global(rho * xr * yr, "sum", na.rm = TRUE)[1, 1])

    list(Ixx = Ixx, Iyy = Iyy, Ixy = Ixy, J = Ixx + Iyy, G = G, rho = rho)
}

#' Concentric-rings J_ref: same construction as gm_depth_index()'s weighted
#' reference (sort by density descending; ring k's outer boundary comes
#' from cumulative cell count, via `.adaptive_density_bins()` in
#' R/utils.R - exact one-ring-per-cell when cheap, a K-bin histogram
#' approximation above that threshold).
#' @param density_sorted numeric vector, density per ring, sorted
#'   descending - from the SAME `.density_raster()`-derived quantity
#'   `.moment_tensor_raster()` used for Ixx/Iyy/Ixy, not
#'   recomputed/renormalised here, to avoid exactly the inconsistency this
#'   file's header documents.
#' @param count numeric vector, number of cells in each ring (default: one
#'   cell per ring)
#' @noRd
.J_ref_raster <- function(cell_area, density_sorted, count = rep(1, length(density_sorted))) {
    S <- cell_area * cumsum(count)
    S_prev <- c(0, S[-length(S)])
    sum(density_sorted * (S^2 - S_prev^2)) / (2 * pi)
}

#' Polar moment-of-inertia compactness/dispersal index for a terra raster -
#' the raster analogue of `shapeindices::moment_of_inertia_index()`. index
#' in `(0, 1]`, `= 1` iff the shape is (almost everywhere) a disk.
#' @inheritParams gm_depth_index
#' @return list(index, J, Ixx, Iyy, Ixy, J_ref, area, centroid)
#' @export
gm_moment_of_inertia_index <- function(rast, weighted = TRUE, n_bins = 1000) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_moment_of_inertia_index")
    valid <- .valid_cells(rast)
    w <- .mass_raster(rast, valid, weighted)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, J = NA_real_, Ixx = NA_real_, Iyy = NA_real_,
                    Ixy = NA_real_, J_ref = NA_real_, area = 0, centroid = NULL))
    }

    core <- .moment_tensor_raster(valid, w, cell_area)
    bins <- .adaptive_density_bins(valid, core$rho, n_bins)
    J_ref <- .J_ref_raster(cell_area, bins$density, bins$count)

    index <- if (core$J > 0) J_ref / core$J else NA_real_
    list(index = index, J = core$J, Ixx = core$Ixx, Iyy = core$Iyy, Ixy = core$Ixy,
         J_ref = J_ref, area = area, centroid = core$G)
}

#' Moment isotropy index for a terra raster - the raster analogue of
#' `shapeindices::moment_isotropy_index()`. index in `(0, 1]`, ratio of
#' the smaller to larger principal moment of the mass inertia tensor -
#' `= 1` iff the mass distribution is rotationally isotropic about its
#' own centroid (any shape with 3-fold or higher rotational symmetry
#' qualifies, not only a disk).
#' @inheritParams gm_depth_index
#' @return list(index, Ixx, Iyy, Ixy, centroid)
#' @export
gm_moment_isotropy_index <- function(rast, weighted = TRUE) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_moment_isotropy_index")
    valid <- .valid_cells(rast)
    w <- .mass_raster(rast, valid, weighted)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, Ixx = NA_real_, Iyy = NA_real_, Ixy = NA_real_, centroid = NULL))
    }

    core <- .moment_tensor_raster(valid, w, cell_area)
    tr <- core$Ixx + core$Iyy
    det <- core$Ixx * core$Iyy - core$Ixy^2
    disc <- sqrt(max(tr^2 - 4 * det, 0))
    lambda_max <- (tr + disc) / 2
    lambda_min <- (tr - disc) / 2

    list(index = if (lambda_max > 0) lambda_min / lambda_max else NA_real_,
         Ixx = core$Ixx, Iyy = core$Iyy, Ixy = core$Ixy, centroid = core$G)
}

#' Directional balance index for a terra raster - the raster analogue of
#' `shapeindices::directional_balance_index()`. index in `[0, 1]`,
#' `= 1 - R` where `R` is the mean resultant length of the mass
#' distribution's own bearing (not distance) from its centroid. `R = 0`
#' (index = 1) means the directional pulls cancel - true for a disk, but
#' also for any shape with 2-fold or higher rotational symmetry (a symmetric
#' dumbbell scores 1 here).
#' @inheritParams gm_depth_index
#' @return list(index, R, centroid)
#' @export
gm_directional_balance_index <- function(rast, weighted = TRUE) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_directional_balance_index")
    valid <- .valid_cells(rast)
    w <- .mass_raster(rast, valid, weighted)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, R = NA_real_, centroid = NULL))
    }

    cell_area <- prod(terra::res(rast))
    cc <- .coord_rasters(valid)
    rho <- .density_raster(valid, w, cell_area)
    G <- .mass_centroid_raster(valid, rho, cc)
    Wr <- as.numeric(terra::global(rho, "sum", na.rm = TRUE)[1, 1])

    xr <- terra::ifel(valid, cc$x - G["x"], NA)
    yr <- terra::ifel(valid, cc$y - G["y"], NA)
    theta <- terra::atan2(yr, xr)

    # R is a MEAN (divided by Wr, rho's own total) - scale-invariant to any
    # constant overall factor in rho, so using the same rho as Ixx/Iyy/Ixy
    # here is for consistency/one-source-of-truth, not because it changes
    # the result versus the raw weight this used before.
    Rx <- as.numeric(terra::global(rho * cos(theta), "sum", na.rm = TRUE)[1, 1]) / Wr
    Ry <- as.numeric(terra::global(rho * sin(theta), "sum", na.rm = TRUE)[1, 1]) / Wr
    R <- sqrt(Rx^2 + Ry^2)

    list(index = 1 - R, R = R, centroid = G)
}
