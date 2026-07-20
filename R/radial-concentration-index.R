## =========================================================================
## Radial concentration index: mean distance to the geometric median vs an
## equal-area circle, for a terra raster - the raster analogue of
## shapeindices::radial_concentration_index().
##
## D1(rho, c) = mean distance from a rho/W-weighted random point to a
## candidate centre c. D1(rho) = min_c D1(rho, c), achieved at the
## geometric median (Fermat-Weber point) - NOT the centroid, which only
## minimises the SQUARED distance (that's gm_moment_of_inertia_index()'s
## job). index = D1_ref/D1(rho), in (0, 1]. See shapeindices' own
## radial-concentration-index.R header for the bathtub-principle proof
## that a disk (unweighted) / concentric rings (weighted) are the exact
## global minimisers - unchanged here, a fact about the density, not
## about triangles vs. raster cells.
##
## ACTUAL VALUE: Weiszfeld's algorithm on the sampled point cloud finds
## the geometric median, then D1 is the mean distance from the sample to
## that point. UNWEIGHTED Weiszfeld here (every sampled point counts
## equally), unlike shapeindices' own weighted `.geometric_median(p, w)` -
## because the density is already baked into WHERE the points were
## sampled from (.sample_valid_points(), weighted or uniform), the same
## importance-sampling argument gm_convexity_index()/gm_span_index() rely on:
## re-weighting an already-density-weighted sample would double-count the
## density.
##
## REFERENCE: D1_ref reuses .annulus_reference_D1_raster() (R/utils.R) -
## the EXACT SAME ring construction gm_depth_index() uses for its own
## weighted reference (there, as `R - D1_ref`; here, directly). No
## quadrature involved (unlike gm_span_index()'s D_ref), so the default
## `n_bins = 1000` (not gm_span_index()'s reduced 100) is appropriate here,
## same reasoning as gm_depth_index()'s own default. ALWAYS this same
## construction, `weighted = FALSE` included: a constant-1 density
## collapses it to `(2/3)*sqrt(area/pi)` exactly (a single degenerate
## ring: r_lo = 0, r_hi = R, mean_r = (2/3)*R, W = 1 - identical algebra to
## `gm_depth_index()`'s file header working through the same collapse for its
## own reference), so there's no separate closed-form path to maintain
## here either.

#' Weiszfeld's algorithm for the (unweighted) geometric median of a point
#' cloud - the point minimising the mean distance to it. Unweighted
#' because the sample itself already reflects the density it was drawn
#' from (see file header) - shapeindices' own version additionally takes
#' a weight vector; not needed here.
#' @param p an Nx2 coordinate matrix
#' @param max_iter maximum iterations
#' @param tol relative convergence tolerance on the objective (mean
#'   distance), not on the centre's position - see shapeindices'
#'   radial-concentration-index.R header for why (a symmetric multi-part
#'   shape's minimiser can be a whole segment, not a single point, so
#'   position can keep moving after the value has converged)
#' @return list(center, D1) - D1 is the converged mean distance
#' @noRd
.geometric_median <- function(p, max_iter = 200, tol = 1e-10) {
    mean_dist <- function(c) mean(sqrt((p[, 1] - c[1])^2 + (p[, 2] - c[2])^2))

    c_t <- colMeans(p)   # start from the centroid
    obj <- mean_dist(c_t)
    for (i in seq_len(max_iter)) {
        d <- pmax(sqrt((p[, 1] - c_t[1])^2 + (p[, 2] - c_t[2])^2), 1e-12)
        c_t <- colSums(p / d) / sum(1 / d)
        obj_new <- mean_dist(c_t)
        if (abs(obj - obj_new) < tol * max(1, obj_new)) { obj <- obj_new; break }
        obj <- obj_new
    }
    list(center = c_t, D1 = obj)
}

#' Radial concentration index for a terra raster - the raster analogue of
#' `shapeindices::radial_concentration_index()`. `D1_ref/D1`, in
#' `(0, 1]`, where `D1` is the mean distance from random interior points
#' to the shape's own geometric median (the point minimising that mean
#' distance - not the centroid) and `D1_ref` is the same quantity for the
#' reference shape (a circle, unweighted; a concentric annulus, weighted) -
#' both provable minimisers of `D1`.
#' @inheritParams gm_depth_index
#' @param size number of points to sample and run Weiszfeld's algorithm
#'   on directly (no pairing, unlike `gm_convexity_index()`/`gm_span_index()`) -
#'   matches `terra::spatSample()`'s own argument name and meaning
#'   directly. Checked against a memory-derived ceiling before running;
#'   hard-stops, not a silent clamp, if exceeded.
#' @param seed optional RNG seed
#' @param n_bins integer, the exact/binned threshold and bin count for the
#'   weighted reference's concentric-rings construction. Default `1000`,
#'   same as `gm_depth_index()`'s - this reference is an exact closed-form
#'   sum, not quadrature, so (unlike `gm_span_index()`) there's no O(rings^2)
#'   cost to guard against.
#' @return `list(index, D1, D1_ref, area, center, n_valid_cells)`. `center`
#'   is the geometric median found by Weiszfeld's algorithm - may be
#'   non-unique for a symmetric multi-part shape (see
#'   `.geometric_median()`'s own comments), in which case it can land
#'   anywhere along the minimising segment, including inside a hole or the
#'   gap between multi-part pieces; the index value itself is unaffected.
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_radial_concentration_index(r, size = 2000, seed = 1)$index
#' @export
gm_radial_concentration_index <- function(rast, weighted = TRUE, size = 3000, seed = NULL, n_bins = 1000) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_radial_concentration_index")
    valid <- .valid_cells(rast)
    w <- .mass_raster(rast, valid, weighted)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, D1 = NA_real_, D1_ref = NA_real_, area = 0,
                    center = NULL, n_valid_cells = 0L))
    }

    .check_mc_size(size, valid, formula = "point", fn_name = "gm_radial_concentration_index")

    if (!is.null(seed)) set.seed(seed)
    pts <- .sample_valid_points(valid, w, size)
    gm <- .geometric_median(pts)

    bins <- .adaptive_density_bins(valid, w, n_bins)
    D1_ref <- .annulus_reference_D1_raster(cell_area, bins$density, bins$count)

    list(index = D1_ref / gm$D1, D1 = gm$D1, D1_ref = D1_ref, area = area,
         center = gm$center, n_valid_cells = as.integer(n_valid))
}
