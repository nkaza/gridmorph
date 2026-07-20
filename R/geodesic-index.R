## =========================================================================
## Geodesic-distance indices: gm_geodesic_span_index() (interior point
## pairs) and gm_geodesic_chord_index() (boundary point pairs, = Angel,
## Parent & Civco (2010)'s "Traversal Index") - the same mean-pairwise-
## distance construction gm_span_index() uses, but the distance between
## two points is the shortest path CONFINED TO the shape, not the
## unconstrained Euclidean straight line. Named "chord", not "traversal",
## to match this package's own already-fixed terminology (a chord is a
## pair on the boundary; span_index()'s own construction uses interior
## points, which is why THAT index isn't called a chord index) - see
## `shapeindices/explorations/NOTES-future-indices.md` item 3, which this
## picks up. That file deferred both indices in the VECTOR package because
## computing a shortest path confined to an arbitrary polygon needs a
## visibility graph (reflex vertices, CGAL, `pathroutr`) - genuinely hard
## computational geometry. NONE of that is needed on a raster.
##
## `terra::gridDist()` (renamed from `gridDistance()` - the old name now
## just errors pointing here) computes exactly the needed distance: for
## every cell, the shortest path to the nearest `target`-valued cell,
## travelling only through non-`NA` cells - a compiled wavefront/cost-
## distance algorithm, already in a dependency this package already uses.
## Verified directly before trusting it (see
## `explorations/NOTES-geodesic-indices.md` for the full record): (1)
## correctly detours around a U-shaped obstacle (2.6x the straight-line
## distance when the direct path is blocked); (2) no directional/
## Manhattan-style bias on open ground (exact match on both diagonal and
## axis-aligned test rays; <=8% bounded angular-quantization error
## elsewhere, the ordinary residual any grid wavefront method has); (3)
## unlike `terra::distance()`/`terra::as.contour()` (`.depth_field()`'s
## and `.contour_perimeter()`'s own padding fixes, R/depth-index.R and
## R/classical-metrics.R), `gridDist()` needs NO edge-padding fix - a
## shape flush against the raster's own true edge gives byte-identical
## results to the same shape with a margin, since it's a pure discrete
## graph computation over cells that actually exist, never needing to
## look past the raster's own extent for anything; (4) converges cleanly
## even on a genuinely serpentine 10-turn test corridor at the default
## `maxiter` - not iteration-starved for any realistic shape.
##
## BOTH INDICES NEED A RESOLUTION-MATCHED REFERENCE - a further, NEW
## finding beyond `vignette("d-resolution-matched-references")`'s own
## four indices, discovered while building these two. `gridDist()`'s own
## angular-quantization error doesn't average away with resolution the
## way `gm_depth_index()`'s cell-centre-vs-edge bias does: it's a property
## of the discrete grid's own connectivity pattern, which repeats
## identically at every scale. Checked directly, not assumed: mean
## geodesic distance between interior point pairs on a rasterized disk,
## compared against `gm_span_index()`'s own closed-form
## `.disk_reference_D()`, stays at ratio ~0.94-0.95 across radius 15/30/60
## (NOT converging to 1 with resolution) - the same non-vanishing-bias
## signature `gm_polsby_popper_index()` had, confirmed for boundary pairs
## too (ratio 1.02-1.04 there). So neither index compares against a pure
## Euclidean-geometry closed form - both build an ACTUAL same-resolution
## reference disk and run it through the identical `gridDist()`-based
## pipeline, exactly the `gm_polsby_popper_index()`/`gm_reock_index()`/
## `gm_detour_index()` pattern (`.reference_disk_raster()`,
## R/utils.R), extended here to a Monte Carlo (not closed-form) reference
## since there's no cheap deterministic geodesic formula to fall back on
## even for a disk.
##
## WEIGHTED VARIANT: span has one, chord does not - a solid disk with
## radial density stays convex regardless of how its own internal density
## varies, so this doesn't change which reference shape is optimal
## (`.reference_disk_raster_weighted()` below reuses the exact same
## cumulative-fraction ring-assignment `.resolution_matched_ref_depth()`
## uses, ranked by radius instead of depth). Boundary cells carry no
## interior mass to weight by, so `gm_geodesic_chord_index()` has no
## `weighted` argument at all, matching the six classical metrics'
## convention rather than the "silently ignored" one.
##
## SAMPLING: K POINTS, NOT K/2 PAIRS - a real efficiency difference from
## every other Monte Carlo index in this package. `gridDist()`'s own
## structure (one call gives distance from ONE source to every other
## cell) makes the usual "size/2 independent pairs" design wasteful here:
## it would need one whole-raster `gridDist()` call PER PAIR. Instead,
## draw K points, run `gridDist()` once per point (K calls), and read off
## the resulting field at the other K-1 sampled points - K(K-1)/2
## pairwise geodesic distances for K distance-transform sweeps. The
## sample mean over all K(K-1)/2 pairs from K i.i.d. draws is still an
## unbiased estimator of E[d(X,Y)] for i.i.d. X,Y (a U-statistic, the same
## structure any all-pairs mean-distance estimator has). Verified fast at
## realistic sizes: K=100 on a 1.44M-cell raster (4,950 pairs) is ~7s,
## scaling roughly linearly in cell count per call.
##
## `n_points`, NOT `size` - A DELIBERATELY DIFFERENT ARGUMENT NAME FROM
## EVERY OTHER MONTE CARLO INDEX IN THIS PACKAGE, and a real design
## reversal worth recording. These two indices were originally built with
## a `size` argument too, matching `gm_span_index()`'s own naming - but
## `gm_span_index(size = N)` draws N points forming N/2 independent
## pairs, while `gm_geodesic_span_index(size = K)` drew K points forming
## ALL K(K-1)/2 pairs among them: the SAME argument name meaning something
## fundamentally different. Once `gm_shape_indices()` was changed to
## include these two in `which = "all"` (see R/shape-indices.R's own
## header), that collision stopped being just a documentation footnote
## and became a real bug: `gm_shape_indices(rast, size = 3000)` - an
## entirely ordinary call, sizing the five PRE-EXISTING Monte Carlo
## indices sensibly - would silently ALSO feed `size = 3000` into these
## two through the shared `...` mechanism, meaning 3000 sequential
## whole-raster `gridDist()` calls each. Verified directly: this measurably
## slowed down an otherwise-ordinary `gm_shape_indices()` call before the
## rename. `n_points` decouples the two families completely - a caller
## tuning `size` for the original five never touches these two by
## accident, and vice versa.
##
## NOISE, HONESTLY: even at `n_points = 200` (verified via repeated
## seeds), mean geodesic distance on a fixed disk carries several percent
## of Monte Carlo noise - genuinely more than `gm_span_index()`'s own
## Euclidean estimate at a comparable sample size, since `gridDist()`'s
## own angular quantization adds variability on top of ordinary point-
## sampling noise. Both the shape's own D and its reference D_ref carry
## this noise independently (drawn from the same `seed`-derived RNG
## stream, but two separate point clouds on two different rasters) -
## `n_points` needs to be noticeably larger here than the Euclidean
## indices' own `size` defaults to get a comparably precise answer, not
## just enough to avoid raw sampling bias.
##
## COST CEILING: `formula = "geodesic"` in `.safe_mc_size_ceiling()`
## (R/utils.R) - genuinely different cost shape from `"point"`/`"line"`
## (K whole-raster `gridDist()` sweeps, each already a chunk-safe,
## terra-native operation on its own - the real resource being protected
## is wall-clock TIME across K sequential calls, not peak memory the way
## the other two formulas are - see that function's own header for the
## full reasoning).
##
## DISCONNECTED SHAPES: detected UP FRONT, before sampling, not just
## after. Geodesic distance is undefined (not just hard to estimate)
## between two points in different connected components of `valid` - a
## genuine limitation, not a computational shortfall (see
## `explorations/NOTES-geodesic-indices.md` for two tempting fixes that
## were checked empirically and rejected: dropping unreachable pairs from
## the mean inflates the index as fragmentation increases - 1.05 to 4.96
## across 1 to 25 equal-area patches in one test - and computing the
## index per patch then averaging discards dispersal information
## entirely, giving IDENTICAL scores whether the same patches sit
## clustered together or scattered across the raster). `.is_connected()`
## checks this with a single `terra::patches(directions = 8)` sweep -
## `directions = 8` because `gridDist()` itself allows diagonal moves
## (verified directly: cost `sqrt(2)` for a diagonal step, confirmed via
## a minimal diagonal-chain test - a `directions = 4` check would
## misclassify a shape connected only through a diagonal pinch point as
## disconnected). Both `gm_geodesic_span_index()` and
## `gm_geodesic_chord_index()` call this BEFORE running the expensive
## K-point `gridDist()` loop, so a disconnected shape returns `NA` with a
## warning at the cost of one cheap raster sweep, not K wasted whole-
## raster sweeps. `.geodesic_pairwise_mean()`'s own per-pair `NaN`
## tracking and `.warn_if_unreachable()` stay in place as a defensive
## backstop - now unreachable in ordinary use, since a confirmed single
## connected component can never produce an unreachable sampled pair, but
## costs nothing to keep.

#' Mean pairwise GEODESIC distance among a set of points, via K
#' `terra::gridDist()` calls (one per point as the sole source, read off
#' at the other K-1 points) - see file header for why this design beats
#' one `gridDist()` call per independent pair.
#'
#' `cells` is DE-DUPLICATED first (by cell number) - unlike
#' `gm_span_index()`'s own size/2-independent-pairs design, EVERY pair
#' among the K points is used here, so a duplicate point would
#' contribute a spurious distance-0 "pair" (itself) that a fresh
#' independent draw would essentially never produce. Achieved K (after
#' de-duplication) is returned alongside D so callers can detect and
#' handle a materially-reduced sample.
#'
#' UNREACHABLE PAIRS: `terra::gridDist()` returns `NaN` (confirmed
#' directly, not assumed) for a cell with no path to the source at all -
#' i.e. the two points sit in different, disconnected components of
#' `valid`. Geodesic distance genuinely has no defined value there (there
#' is no path, not just a long one), unlike Euclidean distance, which is
#' always finite regardless of what lies between the two points -
#' `gm_span_index()` has no analogous failure mode. Rather than let that
#' `NaN` propagate silently into the mean (which it otherwise would, the
#' instant a single sampled pair straddles two components - confirmed
#' directly on a genuinely multi-part raster, not a hypothetical), the
#' count of unreachable pairs actually hit is tracked and returned
#' alongside `D`, so callers can warn with a specific, honest reason
#' instead of surfacing a bare, unexplained `NaN`.
#' @param valid logical terra SpatRaster (TRUE = inside the shape,
#'   traversable)
#' @param cells integer vector of cell numbers, already confirmed to lie
#'   on valid cells
#' @return `list(D, K, n_unreachable_pairs)` - `D` is `NA_real_` if fewer
#'   than 2 unique cells remain after de-duplication, or `NaN` if any
#'   sampled pair turned out to be unreachable (`n_unreachable_pairs > 0`)
#' @noRd
.geodesic_pairwise_mean <- function(valid, cells) {
    cells <- unique(cells)
    K <- length(cells)
    if (K < 2) return(list(D = NA_real_, K = K, n_unreachable_pairs = 0L))

    base <- terra::ifel(valid, 0, NA)
    base_v <- as.vector(terra::values(base))
    total <- 0
    n_unreachable <- 0L
    for (i in seq_len(K - 1)) {
        src_v <- base_v
        src_v[cells[i]] <- 1
        src <- base
        terra::values(src) <- src_v
        d <- terra::gridDist(src, target = 1)
        d_v <- as.vector(terra::values(d))
        pair_d <- d_v[cells[(i + 1):K]]
        n_unreachable <- n_unreachable + sum(is.nan(pair_d))
        total <- total + sum(pair_d)
    }
    list(D = total / (K * (K - 1) / 2), K = K, n_unreachable_pairs = n_unreachable)
}

#' Boundary cells of `valid` - a valid cell adjacent (8-connectivity) to
#' at least one invalid cell OR the raster's own true edge. `fillvalue =
#' 0` in the focal window already correctly treats "beyond the raster's
#' own extent" as invalid for this purpose - verified directly (see file
#' header) that this needs no separate edge-padding fix, unlike
#' `.depth_field()`/`.contour_perimeter()`.
#' @param valid logical terra SpatRaster (TRUE = inside the shape)
#' @return logical terra SpatRaster, TRUE at boundary cells only
#' @noRd
.boundary_cells <- function(valid) {
    interior <- terra::focal(terra::ifel(valid, 1, 0), w = 3, fun = "min", fillvalue = 0)
    valid & (interior == 0)
}

#' Whether `valid`'s TRUE cells form a single connected component, using
#' the SAME connectivity `terra::gridDist()` itself uses - checked
#' directly, not assumed: `gridDist()` allows diagonal moves (cost
#' `sqrt(2)`, confirmed via a minimal diagonal-chain test), i.e.
#' 8-connectivity/Queen's case, so `directions = 4` here would be WRONG -
#' it would flag a shape connected only through a diagonal "pinch point"
#' as disconnected when `gridDist()` itself finds it perfectly reachable.
#'
#' A single `terra::patches()` sweep is `O(valid)` - cheap relative to the
#' K sequential whole-raster `gridDist()` calls `.geodesic_pairwise_mean()`
#' needs. Call this FIRST, before running that expensive loop at all,
#' rather than discovering the same disconnected-shape fact only after
#' paying for it (which is what a pure post-hoc NaN-count check, as in
#' `.geodesic_pairwise_mean()`/`.warn_if_unreachable()`, would do on its
#' own). That post-hoc check stays in place regardless, as a cheap
#' defensive backstop - once this function confirms single-component
#' connectivity under gridDist()'s own rule, every pair within `valid` is
#' reachable by definition of "connected component," so it should never
#' fire again, but costs nothing to keep.
#'
#' Comparing `min == max` of the patch-ID range (rather than counting
#' distinct IDs or taking the max alone) sidesteps `terra::patches()`'s
#' own `allowGaps = TRUE` default, where IDs need not be contiguous
#' (1, 2, 5 is valid) - `min == max` is a robust "at most one patch"
#' check regardless of how IDs happen to be numbered, and only pulls two
#' summary values out of the raster, not full label data.
#' @param valid logical terra SpatRaster (TRUE = inside the shape)
#' @return `TRUE` if `valid` has zero or one connected components (under
#'   8-connectivity), `FALSE` if more than one
#' @noRd
.is_connected <- function(valid) {
    marked <- terra::ifel(valid, 1, NA)
    p <- terra::patches(marked, directions = 8, values = FALSE)
    rng <- terra::global(p, "range", na.rm = TRUE)
    isTRUE(rng[1, 1] == rng[1, 2])
}

#' Warn that the shape itself has more than one connected part, BEFORE any
#' expensive sampling has run (see `.is_connected()`'s own header for why
#' this check comes first) - a specific, honest reason, not a bare `NaN`
#' surfacing later. Distinct from `.warn_if_unreachable()`, which handles
#' the (now purely defensive) case where sampling ran anyway and hit an
#' unreachable pair - this one fires instead of ever starting that loop.
#' @param fn_name caller's own name, for the warning text
#' @return invisible `NULL`
#' @noRd
.warn_disconnected <- function(fn_name) {
    warning(fn_name, "(): the shape has more than one connected part ",
            "(different, disconnected pieces) - mean geodesic distance is ",
            "not defined across disconnected parts, so the index is not ",
            "defined for this shape. Unlike gm_span_index()'s own Euclidean ",
            "distance (always finite regardless of what lies between two ",
            "points), this has no simple fix - a larger `n_points` makes an ",
            "unreachable pair MORE likely to be sampled, not less, for a ",
            "genuinely multi-part shape.")
    invisible(NULL)
}

#' Resolution-matched WEIGHTED reference for `gm_geodesic_span_index()`:
#' an actual rasterized disk at the input's own cell size and area, with
#' cells assigned density by cumulative fraction of RADIUS rank
#' (innermost = highest density) - the same optimal concentric-rings
#' arrangement `gm_span_index()`'s own closed-form reference assumes, but
#' built as a real raster + weight field so it can be sampled through the
#' identical `gridDist()` pipeline the input itself uses. Ranking by
#' radius directly (not via a distance-transform field) is safe here
#' because we build this disk ourselves and already know its own centre
#' exactly - mirrors `.resolution_matched_ref_depth()`'s (R/depth-index.R)
#' cumulative-fraction ring-assignment construction exactly, ranked by
#' radius instead of depth.
#' @param cell_area scalar
#' @param density_sorted numeric vector, density per ring, sorted
#'   descending
#' @param count numeric vector, number of cells in each ring
#' @return `list(valid, weight)`, both terra SpatRasters on the same
#'   (new, disk-sized) grid
#' @noRd
.reference_disk_raster_weighted <- function(cell_area, density_sorted, count = rep(1, length(density_sorted))) {
    n_total <- sum(count)
    ref_valid <- .reference_disk_raster(cell_area, n_total)

    xy <- terra::xyFromCell(ref_valid, seq_len(terra::ncell(ref_valid)))
    valid_v <- as.vector(terra::values(ref_valid))
    ext <- terra::ext(ref_valid)
    cx <- (ext$xmin + ext$xmax) / 2
    cy <- (ext$ymin + ext$ymax) / 2
    r_valid <- sqrt((xy[valid_v, 1] - cx)^2 + (xy[valid_v, 2] - cy)^2)
    n_ref <- length(r_valid)

    cum_frac <- cumsum(count) / n_total
    frac <- (rank(r_valid, ties.method = "first") - 0.5) / n_ref
    ring_idx <- pmin(findInterval(frac, cum_frac) + 1L, length(density_sorted))

    w_v <- rep(NA_real_, length(valid_v))
    w_v[valid_v] <- density_sorted[ring_idx]
    w <- ref_valid
    terra::values(w) <- w_v
    list(valid = ref_valid, weight = w)
}

#' Warn with a specific, honest reason (not a bare `NaN`) when a sampled
#' pair turned out to be unreachable - see `.geodesic_pairwise_mean()`'s
#' own header for why this happens (genuinely disconnected components,
#' not a bug) and converts the resulting `D`/`K`-summary's `NaN` to this
#' package's own `NA_real_` "not defined" convention for consistency with
#' every other undefined-result case elsewhere in this package.
#'
#' Purely a defensive backstop now that both exported functions call
#' `.is_connected()` up front and return via `.warn_disconnected()` before
#' ever sampling - once that check confirms the shape is a single
#' connected component under `gridDist()`'s own connectivity rule, no
#' sampled pair can be unreachable, so this should never actually fire for
#' the `"shape"` side either. Kept anyway rather than removed: cheap, and
#' guards against any mismatch between `.is_connected()`'s assumptions and
#' `gridDist()`'s own behaviour that isn't currently known.
#' @param res a `.geodesic_pairwise_mean()` result
#' @param fn_name caller's own name, for the warning text
#' @param which `"shape"` or `"reference"` - which side hit this, also
#'   for the warning text (the reference disk is always one connected
#'   piece by construction, so this should not normally trigger there,
#'   but is checked anyway rather than assumed)
#' @return `res$D`, with `NaN` converted to `NA_real_`
#' @noRd
.warn_if_unreachable <- function(res, fn_name, which) {
    if (res$n_unreachable_pairs > 0) {
        warning(fn_name, "(): ", res$n_unreachable_pairs, " of the sampled ",
                which, " point pairs were geodesically unreachable from each other ",
                "(different, disconnected parts of the shape) - mean geodesic distance ",
                "is not defined across disconnected parts, so the index is not defined ",
                "for this shape. Unlike gm_span_index()'s own Euclidean distance (always ",
                "finite regardless of what lies between two points), this has no simple ",
                "fix - a larger `n_points` makes an unreachable pair MORE likely to be ",
                "sampled, not less, for a genuinely multi-part shape.")
        return(NA_real_)
    }
    res$D
}

## -- the two indices themselves ------------------------------------------

#' Geodesic span index for a terra raster
#'
#' Mean GEODESIC distance between two random interior points, relative to
#' an equal-area circle - the raster analogue of a shape-confined-path
#' variant of `shapeindices::span_index()` (see file header; deferred in
#' the vector package for lacking cheap shortest-path machinery, not
#' needed here). `index = D_ref/D`, in `(0, 1]`, `= 1` iff the shape is
#' (almost everywhere) a disk. Distinct from `gm_span_index()` (same
#' construction, Euclidean straight-line distance instead) - the two agree
#' exactly on any CONVEX shape (where the straight line between two
#' interior points never leaves it) and diverge only once concavity forces
#' a detour.
#'
#' `n_points` means something DIFFERENT here than `size` does in
#' `gm_span_index()` - a DELIBERATELY DIFFERENT ARGUMENT NAME, not an
#' oversight, precisely to avoid this - see file header. `n_points = K`
#' draws K points and averages ALL K(K-1)/2 pairwise distances among them
#' (not K/2 independent pairs), so a much smaller `n_points` already
#' gives a comparable number of pairs - but each point costs one
#' whole-raster `terra::gridDist()` call, a substantially higher
#' per-point cost than `gm_span_index()`'s own closed-form Euclidean
#' distance. Checked against a memory/time-derived ceiling before running
#' (`formula = "geodesic"`, R/utils.R); hard-stops, not a silent clamp, if
#' exceeded.
#'
#' BOTH `D` and `D_ref` carry real Monte Carlo noise - more than
#' `gm_span_index()`'s own estimate at a comparable sample size, since
#' `terra::gridDist()`'s own angular quantization adds variability on top
#' of ordinary point-sampling noise (verified: several percent spread
#' across seeds even at `n_points = 200` on a test disk). Increase
#' `n_points` for a more precise answer; there is no way to eliminate
#' this noise entirely the way the closed-form Euclidean reference has
#' none.
#' @inheritParams gm_depth_index
#' @param n_points number of points to sample (ALL pairs among them are
#'   used - see above, this is NOT `gm_span_index()`'s own `size`
#'   convention, deliberately a different argument name so the two never
#'   collide through `gm_shape_indices()`'s own shared `...`). Checked
#'   against a memory/time-derived ceiling before running.
#' @param seed optional RNG seed
#' @param n_bins integer, the exact/binned threshold and bin count for the
#'   weighted reference's concentric-rings construction - see
#'   `gm_span_index()`'s own doc for the same parameter's meaning
#' @return `list(index, D, D_ref, area, n_valid_cells)`
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_geodesic_span_index(r, n_points = 25, seed = 1)$index
#' @export
gm_geodesic_span_index <- function(rast, weighted = TRUE, n_points = 40, seed = NULL, n_bins = 100) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_geodesic_span_index")
    valid <- .valid_cells(rast)
    w <- .mass_raster(rast, valid, weighted)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid < 2) {
        if (n_valid == 0) warning("No valid cells; index is not defined.")
        else warning("Only one valid cell; no distinct point pairs to sample - index is not defined.")
        return(list(index = NA_real_, D = NA_real_, D_ref = NA_real_,
                    area = area, n_valid_cells = as.integer(n_valid)))
    }

    if (!.is_connected(valid)) {
        .warn_disconnected("gm_geodesic_span_index")
        return(list(index = NA_real_, D = NA_real_, D_ref = NA_real_,
                    area = area, n_valid_cells = as.integer(n_valid)))
    }

    .check_mc_size(n_points, valid, formula = "geodesic", fn_name = "gm_geodesic_span_index", arg_name = "n_points")
    if (!is.null(seed)) set.seed(seed)

    pts <- .sample_valid_points(valid, w, max(2L, n_points))
    cells <- terra::cellFromXY(valid, pts)
    actual <- .geodesic_pairwise_mean(valid, cells)

    bins <- .adaptive_density_bins(valid, w, n_bins)
    ref <- .reference_disk_raster_weighted(cell_area, bins$density, bins$count)
    ref_pts <- .sample_valid_points(ref$valid, ref$weight, max(2L, n_points))
    ref_cells <- terra::cellFromXY(ref$valid, ref_pts)
    reference <- .geodesic_pairwise_mean(ref$valid, ref_cells)

    D <- .warn_if_unreachable(actual, "gm_geodesic_span_index", "shape")
    D_ref <- .warn_if_unreachable(reference, "gm_geodesic_span_index", "reference")
    index <- if (is.na(D) || is.na(D_ref)) NA_real_ else D_ref / D

    list(index = index, D = D, D_ref = D_ref, area = area, n_valid_cells = as.integer(n_valid))
}

#' Geodesic chord index for a terra raster
#'
#' Mean GEODESIC distance between two random points on the shape's own
#' BOUNDARY, relative to an equal-area circle - the raster analogue of
#' Angel, Parent & Civco (2010)'s "Traversal Index" (see file header for
#' the naming choice and why it's called "chord" here, not "traversal").
#' `index = D_ref/D`, in `(0, 1]`, `= 1` iff the shape is (almost
#' everywhere) a disk.
#'
#' No `weighted` argument - boundary cells carry no interior mass to
#' weight by (see file header), matching `gm_hull_ratio_index()`/
#' `gm_polsby_popper_index()`/etc.'s convention of omitting the argument
#' entirely rather than silently ignoring it.
#'
#' `n_points` points are drawn WITHOUT replacement from the finite set of
#' boundary cells (capped at however many exist, with a warning if
#' `n_points` had to be reduced) - unlike `gm_geodesic_span_index()`'s own
#' with-replacement interior sampling, there is no density to weight by
#' here, so plain uniform sampling over a known finite population is both
#' simpler and avoids any duplicate-point risk by construction. A
#' deliberately different argument name from `gm_span_index()`'s own
#' `size` - see `gm_geodesic_span_index()`'s own doc for why.
#' @inheritParams gm_hull_ratio_index
#' @param n_points number of boundary points to sample (ALL pairs among
#'   them are used - see `gm_geodesic_span_index()`'s own doc for why
#'   this argument is named differently from `gm_span_index()`'s `size`).
#'   Checked against a memory/time-derived ceiling before running.
#' @param seed optional RNG seed
#' @return `list(index, D, D_ref, area, n_valid_cells, n_boundary_cells)`
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_geodesic_chord_index(r, n_points = 25, seed = 1)$index
#' @export
gm_geodesic_chord_index <- function(rast, n_points = 40, seed = NULL) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_geodesic_chord_index")
    valid <- .valid_cells(rast)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, D = NA_real_, D_ref = NA_real_,
                    area = 0, n_valid_cells = 0L, n_boundary_cells = 0L))
    }

    bnd <- .boundary_cells(valid)
    bnd_cells <- which(as.vector(terra::values(bnd)))
    n_bnd <- length(bnd_cells)
    if (n_bnd < 2) {
        warning("Fewer than two boundary cells; index is not defined.")
        return(list(index = NA_real_, D = NA_real_, D_ref = NA_real_, area = area,
                    n_valid_cells = as.integer(n_valid), n_boundary_cells = n_bnd))
    }

    if (!.is_connected(valid)) {
        .warn_disconnected("gm_geodesic_chord_index")
        return(list(index = NA_real_, D = NA_real_, D_ref = NA_real_, area = area,
                    n_valid_cells = as.integer(n_valid), n_boundary_cells = n_bnd))
    }

    .check_mc_size(n_points, valid, formula = "geodesic", fn_name = "gm_geodesic_chord_index", arg_name = "n_points")
    if (n_points > n_bnd) {
        warning("`n_points` (", n_points, ") exceeds the number of boundary cells (", n_bnd,
                "); sampling all ", n_bnd, " instead.")
    }
    k <- min(n_points, n_bnd)
    if (!is.null(seed)) set.seed(seed)

    pick <- sample(bnd_cells, k)
    actual <- .geodesic_pairwise_mean(valid, pick)

    ref_valid <- .reference_disk_raster(cell_area, n_valid)
    ref_bnd <- which(as.vector(terra::values(.boundary_cells(ref_valid))))
    ref_pick <- sample(ref_bnd, min(k, length(ref_bnd)))
    reference <- .geodesic_pairwise_mean(ref_valid, ref_pick)

    D <- .warn_if_unreachable(actual, "gm_geodesic_chord_index", "shape")
    D_ref <- .warn_if_unreachable(reference, "gm_geodesic_chord_index", "reference")
    index <- if (is.na(D) || is.na(D_ref)) NA_real_ else D_ref / D

    list(index = index, D = D, D_ref = D_ref,
         area = area, n_valid_cells = as.integer(n_valid), n_boundary_cells = n_bnd)
}
