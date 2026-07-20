## =========================================================================
## Convexity / dispersal index for a terra raster - the raster analogue of
## shapeindices::convexity_index().
##
## Definition (unchanged from the vector package): CI = 1 - E[f(X,Y)],
## where X, Y are random points drawn from inside the shape and f(X,Y) is
## the fraction of segment XY lying outside the shape. CI = 1 iff the
## shape is convex.
##
## MONTE CARLO ONLY - no deterministic mode, unlike shapeindices'
## deterministic = TRUE/FALSE split. shapeindices' deterministic = TRUE is
## exhaustive all-pairs over CDT triangles, practical there because a real
## mesh is tens to a few hundred triangles. The raster analogue of
## "triangle" is "valid cell" - even a modest 1000-cell raster is ~500K
## pairs, and a 100K-cell raster (unremarkable for a real raster) is ~5e9 -
## exhaustive all-cell-pairs is never practical at realistic raster sizes
## (see the package's own build plan for the full reasoning). Monte Carlo
## is the only mode here, not a fallback.
##
## LINE-OUTSIDE FRACTION: rather than a vector line-clipping computation
## (GEOS calls, boundary-edge intersection), each sampled line is
## discretised directly at the raster's own resolution and looked up
## against `valid` - see .frac_outside_line() in R/utils.R. WEIGHTED
## sampling (denser draws where `weight` is higher, via
## .sample_valid_points()) already produces a correctly weighted
## expectation via a plain mean over sampled lines - no further per-line
## weighting needed, the same importance-sampling argument shapeindices'
## own Monte Carlo mode relies on.

#' Convexity/dispersal index of a terra raster - the raster analogue of
#' `shapeindices::convexity_index()`. 1 minus the expected fraction of a
#' random interior line lying outside the shape. index in `[0, 1]`, `1` =
#' convex; lower means more concave and/or more spatially dispersed.
#' Handles holes and multi-part shapes with no special-casing: a line
#' between two disjoint parts is mostly "outside", exactly like a vector
#' MULTIPOLYGON's own line-clipping.
#' @inheritParams gm_depth_index
#' @param size number of points to sample (paired up consecutively into
#'   `size %/% 2` lines) - matches `terra::spatSample()`'s own argument
#'   name and meaning directly (a raw point count, not a line count).
#'   Checked against a memory-derived ceiling before running (see
#'   `.safe_mc_size_ceiling()` in `R/utils.R`); hard-stops, not a silent
#'   clamp, if exceeded.
#' @param seed optional RNG seed
#' @return `list(index, n_lines, area, n_valid_cells)`
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_convexity_index(r, size = 500, seed = 1)$index
#' @export
gm_convexity_index <- function(rast, weighted = TRUE, size = 3000, seed = NULL) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_convexity_index")
    valid <- .valid_cells(rast)
    w <- .mass_raster(rast, valid, weighted)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, n_lines = 0L, area = 0, n_valid_cells = 0L))
    }
    if (n_valid == 1) {
        return(list(index = 1, n_lines = 0L, area = area, n_valid_cells = 1L))  # vacuously convex
    }

    .check_mc_size(size, valid, formula = "line", fn_name = "gm_convexity_index")

    n_lines <- max(1L, size %/% 2L)
    if (!is.null(seed)) set.seed(seed)
    pts <- .sample_valid_points(valid, w, 2L * n_lines)
    x1 <- pts[seq(1, 2 * n_lines, 2), , drop = FALSE]
    x2 <- pts[seq(2, 2 * n_lines, 2), , drop = FALSE]

    frac_outside <- .frac_outside_line(valid, x1, x2)
    index <- 1 - mean(frac_outside)

    list(index = index, n_lines = n_lines, area = area, n_valid_cells = as.integer(n_valid))
}
