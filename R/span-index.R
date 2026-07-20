## =========================================================================
## Span index: mean pairwise interior distance vs. an equal-area circle,
## for a terra raster - the raster analogue of shapeindices::span_index().
##
## D(rho) = mean distance between two independent points drawn from
## rho/W. index = D_ref/D, in (0, 1], = 1 iff rho already matches its
## reference shape (a disk, unweighted; concentric rings sorted
## densest-to-centre, weighted) - see shapeindices' own span-index.R
## header for the full Riesz-rearrangement-inequality proof, unchanged
## here since it's a fact about the density itself, not about triangles
## vs. raster cells.
##
## MONTE CARLO ONLY, same reasoning as gm_convexity_index()'s own header:
## shapeindices' deterministic = TRUE mode is exhaustive over CDT
## triangle-pairs (tens to a few hundred), impractical over raster
## cell-pairs (a modest 1000-cell raster is already ~500K pairs). No
## same-triangle "self-term" issue to port either - that complication is
## specific to shapeindices' deterministic quadrature-MESH mode
## (.mesh_span_index()), which isn't built here at all: a raster's own
## Monte Carlo point pairs (sampled directly at cell centres) already ARE
## the estimator of D, the same way shapeindices' own
## .random_pair_span_index() needs no self-term either.
##
## ACTUAL VALUE: mean Euclidean distance over `size %/% 2` sampled point
## pairs - see .sample_valid_points() in R/utils.R for the (weighted or
## uniform) sampling itself.
##
## WEIGHTED REFERENCE: no closed form exists ring-to-ring (the |x-y|
## kernel has no polynomial antiderivative in r the way
## gm_moment_of_inertia_index()'s/gm_depth_index()'s r^2 or r do), but averaging
## over the angle between two points at radii r1, r2 does - the complete
## elliptic integral of the second kind, .ellip_E()/.mean_chord() below,
## ported VERBATIM from shapeindices::span-index.R (pure math, no sf/
## triangulation dependency to strip out; identical for a raster's
## concentric rings as for a polygon's).
##
## RING COUNT: reuses .adaptive_density_bins() (R/utils.R) for the same
## exact/binned adaptive construction gm_depth_index()/gm_moment_of_inertia_index()
## use - but with a MUCH SMALLER default n_bins (100, not 1000) here,
## because THIS reference's own cost is O(rings^2 * gl_order^2) (a
## quadrature-node pair sum, not the O(rings) cumulative-sum the other two
## indices' references need) - directly mirrors shapeindices'
## .annulus_reference_D()'s own choice of `max_rings = 100` for the exact
## same reason (one ring per triangle instead crashes past 18GB on a
## real-world mesh - see that function's own comments).
## n_bins = 100 with gl_order = 8 (fixed, not user-facing, matching
## shapeindices' own choice) means at most 800 quadrature nodes and a
## 800x800 pairwise sum - a few MB, not the ~1GB an n_bins = 1000 version
## of this SAME reference would need.

## -- special functions and quadrature, ported verbatim from shapeindices --

#' Complete elliptic integral of the second kind via the AGM algorithm.
#' @param k numeric vector in `[0, 1]`
#' @return E(k), same length as k
#' @noRd
.ellip_E <- function(k) {
    k <- pmin(pmax(k, 0), 1)
    out <- numeric(length(k))
    near1 <- k > 1 - 1e-9
    out[near1] <- 1   # E(1) = 1 exactly; AGM below is 0/Inf indeterminate there
    idx <- which(!near1)
    if (length(idx) == 0) return(out)

    a <- rep(1, length(idx)); b <- sqrt(1 - k[idx]^2); c <- k[idx]
    sum_term <- 0.5 * c^2
    pow2 <- 0.5
    for (i in seq_len(30)) {
        a_new <- (a + b) / 2
        b_new <- sqrt(a * b)
        c <- (a - b) / 2
        a <- a_new; b <- b_new
        pow2 <- pow2 * 2
        sum_term <- sum_term + pow2 * c^2
        if (max(c) < 1e-16) break
    }
    out[idx] <- (pi / (2 * a)) * (1 - sum_term)
    out
}

#' Mean distance between 2 points at radii r1, r2 with uniform angle
#' between them - the angular part of the annulus reference integral.
#' @param r1,r2 numeric vectors (recycled against each other)
#' @return same length as the longer of r1, r2
#' @noRd
.mean_chord <- function(r1, r2) {
    s <- r1 + r2
    k <- ifelse(s > 0, pmin(2 * sqrt(r1 * r2) / s, 1), 0)
    (2 / pi) * s * .ellip_E(k)
}

#' m-point Gauss-Legendre nodes/weights on `[-1, 1]` via the Golub-Welsch
#' eigendecomposition of the Jacobi matrix.
#' @param m number of nodes
#' @return list(x, w)
#' @noRd
.gauss_legendre <- function(m) {
    if (m == 1) return(list(x = 0, w = 2))
    k <- seq_len(m - 1)
    beta <- k / sqrt(4 * k^2 - 1)
    J <- matrix(0, m, m)
    J[cbind(seq_len(m - 1), 2:m)] <- beta
    J[cbind(2:m, seq_len(m - 1))] <- beta
    eig <- eigen(J, symmetric = TRUE)
    ord <- order(eig$values)
    list(x = eig$values[ord], w = 2 * eig$vectors[1, ord]^2)
}

#' Gauss-Legendre nodes/weights transformed to radii in `[r_lo, r_hi]`,
#' weighted by that ring's own radial density f(r) = 2r/(r_hi^2 - r_lo^2)
#' (the radius marginal of a uniform annulus).
#' @param r_lo,r_hi ring bounds
#' @param gl a .gauss_legendre() result, reused across rings
#' @return list(r, p) - p sums to 1
#' @noRd
.radial_nodes <- function(r_lo, r_hi, gl) {
    half <- (r_hi - r_lo) / 2
    mid  <- (r_hi + r_lo) / 2
    r    <- half * gl$x + mid
    dens <- 2 * r / (r_hi^2 - r_lo^2)
    p    <- gl$w * half * dens
    list(r = r, p = p / sum(p))
}

## -- reference D -------------------------------------------------------

#' @noRd
.disk_reference_D <- function(area) {
    (128 / (45 * pi)) * sqrt(area / pi)
}

#' Weighted reference: rings from .adaptive_density_bins() (exact per-cell
#' below n_bins valid cells, a K-bin histogram above), then average
#' .mean_chord() over every pair of rings' own radial marginals -
#' flattened into one radial "point cloud" so it's a single vectorised
#' sum. See file header for why `n_bins` defaults much lower here than in
#' gm_depth_index()/gm_moment_of_inertia_index() (O(rings^2), not O(rings)).
#' @param cell_area scalar, the (identical) area of one raster cell
#' @param density_sorted numeric vector, density per ring, sorted
#'   descending
#' @param count numeric vector, number of cells in each ring
#' @param gl_order Gauss-Legendre nodes per ring (fixed, not user-facing,
#'   matching shapeindices' own choice)
#' @return D_ref
#' @noRd
.annulus_reference_D_raster <- function(cell_area, density_sorted, count = rep(1, length(density_sorted)), gl_order = 8) {
    S    <- cell_area * cumsum(count)
    r_hi <- sqrt(S / pi)
    r_lo <- c(0, r_hi[-length(r_hi)])
    mass <- density_sorted * count
    W    <- mass / sum(mass)

    gl    <- .gauss_legendre(gl_order)
    nodes <- Map(.radial_nodes, r_lo, r_hi, MoreArgs = list(gl = gl))
    r_all <- unlist(lapply(nodes, `[[`, "r"))
    p_all <- unlist(Map(function(nd, w_i) nd$p * w_i, nodes, W))

    sum(outer(p_all, p_all) * outer(r_all, r_all, .mean_chord))
}

## -- the span index itself ----------------------------------------------

#' Span index for a terra raster - the raster analogue of
#' `shapeindices::span_index()`. `D_ref/D`, in `(0, 1]`, where `D` is the
#' mean distance between two random interior points and `D_ref` is the
#' same quantity for the reference shape (a circle, unweighted; a
#' concentric annulus, weighted) - both provable minimisers of `D`.
#' Distinct from `gm_moment_of_inertia_index()`, which a squared-distance
#' version of this index would just collapse to.
#' @inheritParams gm_depth_index
#' @param size number of points to sample (paired up consecutively into
#'   `size %/% 2` point pairs) - matches `terra::spatSample()`'s own
#'   argument name and meaning directly. Checked against a memory-derived
#'   ceiling before running; hard-stops, not a silent clamp, if exceeded.
#' @param seed optional RNG seed
#' @param n_bins integer, the exact/binned threshold and bin count for the
#'   weighted reference's concentric-rings construction. Default `100`,
#'   much lower than `gm_depth_index()`'s/`gm_moment_of_inertia_index()`'s
#'   `1000` - see file header for why (this reference's own cost is
#'   quadratic in ring count, not linear).
#' @return `list(index, D, D_ref, area, n_valid_cells)`
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_span_index(r, size = 2000, seed = 1)$index
#' @export
gm_span_index <- function(rast, weighted = TRUE, size = 3000, seed = NULL, n_bins = 100) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_span_index")
    valid <- .valid_cells(rast)
    w <- .mass_raster(rast, valid, weighted)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, D = NA_real_, D_ref = NA_real_, area = 0, n_valid_cells = 0L))
    }

    .check_mc_size(size, valid, formula = "point", fn_name = "gm_span_index")

    n_pairs <- max(1L, size %/% 2L)
    if (!is.null(seed)) set.seed(seed)
    pts <- .sample_valid_points(valid, w, 2L * n_pairs)
    x1 <- pts[seq(1, 2 * n_pairs, 2), , drop = FALSE]
    x2 <- pts[seq(2, 2 * n_pairs, 2), , drop = FALSE]

    D <- mean(sqrt((x1[, 1] - x2[, 1])^2 + (x1[, 2] - x2[, 2])^2))

    bins <- .adaptive_density_bins(valid, w, n_bins)
    D_ref <- if (length(bins$density) == 1) {
        # a spatially-constant density (includes weighted = FALSE)
        # collapses to one ring covering the whole shape - that IS the
        # unweighted disk case algebraically, so use the closed form
        # directly rather than a single ring's own 8-point self-pairing
        # quadrature, which (verified) is measurably less accurate than
        # the exact formula for exactly this degenerate case: an 8-node
        # rule integrating the WHOLE [0, R] range against itself is far
        # coarser than either the closed form or a many-ring subdivision
        # of the same uniform density would be
        .disk_reference_D(area)
    } else {
        .annulus_reference_D_raster(cell_area, bins$density, bins$count)
    }

    list(index = D_ref / D, D = D, D_ref = D_ref, area = area, n_valid_cells = as.integer(n_valid))
}
