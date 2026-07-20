## =========================================================================
## Six classical, purely boundary/hull/extent-based compactness scores,
## ported from shapeindices::classical-metrics.R - the raster analogue.
## Unrelated to (and much cheaper than) the sampling/moment-based indices
## elsewhere in this package.
##
##   gm_hull_ratio_index         = area / area(convex hull)
##   gm_polsby_popper_index      = 4*pi*area / perimeter^2
##   gm_width_length_ratio_index = min(bbox width, height) / max(...)
##   gm_reock_index               = area / area(minimum bounding circle)
##   gm_detour_index               = perimeter(equal-area circle) / perimeter(convex hull)
##   gm_exchange_index             = area(shape INTERSECT equal-area circle at centroid) / area
##
## NO `weighted` ARGUMENT - unlike every other index in this package,
## these six are structurally unweighted, matching shapeindices' own
## explicit, reasoned decision (see its own R/classical-metrics.R header):
## area, perimeter, convex hull, and minimum bounding circle are all
## properties of the shape's own EXTENT/BOUNDARY, unaffected by any
## redistribution of mass within it - there's nothing for a "weighted"
## version to act on. Substituting a weighted area for
## gm_polsby_popper_index()/gm_reock_index() specifically was considered and
## rejected there for a precise reason: their (0, 1] boundedness comes
## from a theorem tying BOTH sides of the ratio together (isoperimetric
## inequality; MBC containment) - swap one side for an arbitrarily-scaled
## weight and nothing keeps the ratio bounded any more. Same reasoning
## carries over unchanged here, since it's a fact about the theorems, not
## about triangles vs. raster cells.
##
## PIXEL-COUNTING AREA: area = n_valid_cells * cell_area throughout, the
## same convention every other index in this package already uses.
##
## PERIMETER CONVENTION - decided empirically, not by default: naive
## cell-edge/stair-step boundary tracing (terra::as.polygons()'s own
## perimeter) badly OVERESTIMATES the true perimeter for any non-axis-
## aligned boundary - verified directly on two test shapes at matched
## scale: +46% on a 45-degree diamond, +32% on a circle.
## terra::as.contour() (marching squares at the 0.5 level, the 0/1
## midpoint) is substantially closer in both cases (+3.3% and +7.9%
## respectively) - still a real bias (marching squares is itself an
## approximation, not exact), but the clearly better raster-native
## convention, so it's what gm_polsby_popper_index() uses for the shape's
## OWN (non-convex, possibly multi-part) perimeter. gm_detour_index()'s
## hull perimeter has no such issue at all - terra::convHull()'s own
## output is already a proper polygon with straight edges, not a raster-
## traced boundary, so terra::perim() on it is exact regardless of
## orientation.
##
## RASTER-EDGE PADDING: `.contour_perimeter()` (below) pads before calling
## terra::as.contour() - a real, verified bug otherwise, the same
## underlying problem as R/morphology.R's own focal() fillvalue fix and
## R/depth-index.R's own terra::distance() padding fix, just a third
## terra function: marching squares needs cell values on both sides of a
## boundary to interpolate it, and has nothing beyond the raster's own
## true edge to interpolate against there, so a shape's own edge-touching
## side goes untraced. See .contour_perimeter()'s own comment for the
## full reproduction (perimeter roughly HALVED for a square flush against
## the raster's edge, before the fix).
##
## MULTI-PART SHAPES: terra::as.polygons(dissolve = TRUE) combines every
## valid cell (however many disjoint groups) into ONE (possibly multi-
## part) polygon before any hull/perimeter computation - convHull() on
## that spans across every part, matching shapeindices' own "hull over
## the union of all parts" convention with no special-casing needed.
##
## RESOLUTION-MATCHED REFERENCES, HARDWIRED (no opt-out) - three of these
## six scores compare a raster-MEASURED quantity against a pure ANALYTICAL
## formula (a perfect circle's `1`, or its exact perimeter/area), the same
## disease `gm_depth_index()` had (see that file's own header): a
## rasterized disk should score `1` on `gm_polsby_popper_index()`/
## `gm_reock_index()`/`gm_detour_index()`, but scored `0.880`/`0.953`/
## `0.982` respectively at a representative resolution, comparing raw
## against the pure formula. Fixed the same way depth was: build an ACTUAL
## rasterized disk at the input's own cell size and area, compute ITS OWN
## raw score through the IDENTICAL pipeline, divide by that instead of by
## the analytical constant.
##
## `gm_polsby_popper_index()` NEEDS PHASE-AVERAGING; `gm_reock_index()`/
## `gm_detour_index()` DON'T - checked, not assumed, before deciding this.
## `.contour_perimeter()`'s marching-squares perimeter, verified against
## 25 sub-cell offsets per radius: (1) its PHASE sensitivity is large (up
## to several percent at a single realization, vs. `gm_depth_index()`'s
## own ~0.1%), because marching squares interpolates a 0.5 crossing
## against a raster that was never smoothly sampled in the first place -
## there's no genuine sub-cell signal there to interpolate, unlike
## `.depth_field()`'s distance values; (2) its MEAN (properly phase-
## averaged) does NOT converge toward `1` as resolution increases either -
## it plateaus around `0.90` from radius 10 to radius 300 - a fundamentally
## different, non-vanishing bias, not depth's `O(cell_size)` one. A SINGLE
## reference-disk realization is therefore not a reliable reference for
## this one (verified: it left the corrected index drifting as far as
## `1.04` at some radii); a 25-point sub-cell offset grid was verified
## sufficient (corrected index landed within `+-0.006` of `1` across
## radius 30 to 300, using a properly-averaged reference at every one).
## `gm_reock_index()`'s area (pixel count) and `gm_detour_index()`'s hull
## perimeter (an EXACT polygon, no marching-squares reconstruction
## involved at all - see this file's own `gm_detour_index()` roxygen) have
## no such interpolation-noise source, so both behave like
## `gm_depth_index()` did: low phase sensitivity (sd `<0.004` at every
## radius checked) and a reference that converges smoothly toward `1`
## with resolution - a single reference-disk realization is enough for
## both, verified directly (`0.997-1.000` across radius 15 to 120).
##
## `gm_hull_ratio_index()`/`gm_exchange_index()` are DELIBERATELY NOT
## included in this treatment. Both are self-referential (shape vs. its
## OWN convex hull; shape vs. itself intersected with a circle at its own
## centroid) - there is no separate analytical formula being compared
## against in the first place, so there's nothing for a same-resolution
## reference disk to correct. `gm_hull_ratio_index()`'s own sub-1 score
## for a disk (`0.968`) reflects a REAL property of the rasterized shape -
## a pixelated disk's stair-step boundary genuinely isn't perfectly convex
## at the pixel level, the same way a real coastline genuinely isn't a
## smooth curve - not a measurement artifact to correct away.
##
## gm_exchange_index()'s KNOWN PATHOLOGY carries over unchanged from
## shapeindices: the reference circle is centred at the OVERALL
## (unweighted) centroid, which for well-separated multi-part shapes can
## sit in the empty gap between parts - once the circle is smaller than
## the distance to the nearest part, the index is exactly 0. Included
## despite this, matching shapeindices' own precedent of shipping a real,
## useful, documented-blind-spot metric rather than omitting it.

## -- minimum enclosing circle, for gm_reock_index() - ported VERBATIM from
## shapeindices::classical-metrics.R (pure coordinate-geometry math, no
## sf/triangulation dependency to strip out) -------------------------------

#' @noRd
.circle_from_2 <- function(p1, p2) {
    list(center = (p1 + p2) / 2, r = sqrt(sum((p1 - p2)^2)) / 2)
}

#' Circumcircle of 3 points, falling back to the widest pairwise circle if
#' they're (near-)collinear.
#' @noRd
.circle_from_3 <- function(p1, p2, p3) {
    ax <- p1[1]; ay <- p1[2]; bx <- p2[1]; by <- p2[2]; cx <- p3[1]; cy <- p3[2]
    d <- 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if (abs(d) < 1e-9) {
        cands <- list(.circle_from_2(p1, p2), .circle_from_2(p2, p3), .circle_from_2(p1, p3))
        return(cands[[which.max(vapply(cands, `[[`, numeric(1), "r"))]])
    }
    ux <- ((ax^2 + ay^2) * (by - cy) + (bx^2 + by^2) * (cy - ay) + (cx^2 + cy^2) * (ay - by)) / d
    uy <- ((ax^2 + ay^2) * (cx - bx) + (bx^2 + by^2) * (ax - cx) + (cx^2 + cy^2) * (bx - ax)) / d
    center <- c(ux, uy)
    list(center = center, r = sqrt(sum((p1 - center)^2)))
}

#' @noRd
.point_in_circle <- function(p, circle, eps = 1e-7) {
    sqrt(sum((p - circle$center)^2)) <= circle$r + eps
}

#' Minimum enclosing circle of a point set - the classic iterative
#' Welzl-style incremental algorithm, without the usual random shuffle (so
#' it never touches the caller's RNG state); worst case O(n^3), fine for
#' the vertex counts a convex hull actually has.
#' @param pts an Nx2 coordinate matrix
#' @return list(center, r)
#' @noRd
.min_enclosing_circle <- function(pts) {
    n <- nrow(pts)
    if (n == 1) return(list(center = pts[1, ], r = 0))

    circle <- .circle_from_2(pts[1, ], pts[2, ])
    for (i in seq_len(n)) {
        if (.point_in_circle(pts[i, ], circle)) next
        circle <- list(center = pts[i, ], r = 0)
        for (j in seq_len(i - 1)) {
            if (.point_in_circle(pts[j, ], circle)) next
            circle <- .circle_from_2(pts[i, ], pts[j, ])
            for (k in seq_len(j - 1)) {
                if (.point_in_circle(pts[k, ], circle)) next
                circle <- .circle_from_3(pts[i, ], pts[j, ], pts[k, ])
            }
        }
    }
    circle
}

## -- shared raster-native geometry helpers ---------------------------------

#' The valid-cell footprint as ONE (possibly multi-part) polygon, dissolved
#' across every valid cell - matches shapeindices' own union-before-hull
#' convention for multi-part shapes, achieved here for free by
#' `terra::as.polygons(dissolve = TRUE)` rather than an explicit union step.
#' @param valid logical terra SpatRaster (TRUE = inside the shape)
#' @return a terra SpatVector, one row, polygons (possibly multi-part)
#' @noRd
.valid_polygons <- function(valid) {
    p <- terra::as.polygons(terra::ifel(valid, 1, 0), dissolve = TRUE)
    p[p[[1]][[1]] == 1, ]
}

#' Raw shape perimeter via marching squares (`terra::as.contour()` at the
#' 0/1 midpoint) - see this file's own header for why, verified
#' empirically, not the naive cell-edge/stair-step alternative.
#'
#' PADS the 0/1 raster by one cell (`terra::extend(..., fill = 0)`) before
#' contouring - a real, verified bug otherwise, the `as.contour()`
#' analogue of `R/morphology.R`'s own `focal()` `fillvalue` fix and
#' `gm_depth_index()`'s own `terra::distance()` padding fix (three
#' different terra functions, the same underlying problem: marching
#' squares needs cell values on both sides of a boundary to interpolate
#' it, and has nothing beyond the raster's own true edge to interpolate
#' against there). Verified directly: a 21-cell-wide square flush against
#' the raster's top edge gave perimeter `40`, roughly HALF the correct
#' `82.8` (`terra::as.contour()` couldn't trace that side at all without
#' padding) - the padded fix gives EXACTLY the same perimeter as the
#' identical square positioned away from the edge with a margin.
#' @param valid logical terra SpatRaster (TRUE = inside the shape)
#' @return numeric, total perimeter length in `rast`'s own CRS units
#' @noRd
.contour_perimeter <- function(valid) {
    padded <- terra::extend(terra::ifel(valid, 1, 0), 1, fill = 0)
    cont <- terra::as.contour(padded, levels = 0.5)
    # muffled: redundant with .check_planar_crs()'s own, clearer warning
    # for a missing CRS, already issued once at the top of the caller
    sum(.muffle_crs_warnings(terra::perim(cont)))
}

## -- resolution-matched references, see this file's own header ------------

#' `gm_polsby_popper_index()`'s reference: the raw `4*pi*area/perimeter^2`
#' score of an actual rasterized disk at the input's own cell size and
#' area, AVERAGED over a grid of sub-cell centre offsets - see file header
#' for why a single realization isn't good enough here specifically (large
#' phase sensitivity, non-vanishing mean bias) where it's fine for
#' `.disk_raw_detour()`/`.disk_raw_reock()` below.
#'
#' A diagonal offset near the grid's own corners can miss every cell
#' centre entirely when `n_valid` is very small (e.g. `n_valid = 1`: the
#' reference disk's radius is `~0.56` cells, less than the worst-case
#' `~0.64`-cell diagonal offset) - `terra::as.contour()` has nothing to
#' trace on an all-empty raster. Those phases are dropped from the
#' average rather than propagating a crash; at least one phase always
#' succeeds; averaging over fewer, less-degenerate phases is still a
#' better reference than none.
#' @param cell_area scalar
#' @param n_valid target valid-cell count for the reference disk
#' @param n_phase_side sub-cell offsets sampled per axis (25 total at the
#'   default 5 - verified sufficient, see file header)
#' @return mean raw polsby-popper score across the phase grid
#' @noRd
.phase_averaged_ref_pp <- function(cell_area, n_valid, n_phase_side = 5) {
    offsets <- seq(-0.45, 0.45, length.out = n_phase_side)
    phases <- expand.grid(x = offsets, y = offsets)
    scores <- vapply(seq_len(nrow(phases)), function(i) {
        ref_valid <- .reference_disk_raster(cell_area, n_valid, phase = c(phases$x[i], phases$y[i]))
        ref_n <- as.numeric(terra::global(ref_valid, "sum", na.rm = TRUE)[1, 1])
        if (ref_n == 0) return(NA_real_)
        ref_area <- ref_n * cell_area
        ref_perim <- .contour_perimeter(ref_valid)
        if (ref_perim == 0) return(NA_real_)
        (4 * pi * ref_area) / ref_perim^2
    }, numeric(1))
    mean(scores, na.rm = TRUE)
}

#' `gm_detour_index()`'s reference: the raw `circle_perimeter/hull_perimeter`
#' score of an actual rasterized disk at the input's own cell size and
#' area - a single realization, unlike `.phase_averaged_ref_pp()` above
#' (verified low phase sensitivity, see file header).
#' @param cell_area scalar
#' @param n_valid target valid-cell count for the reference disk
#' @return raw detour score of the reference disk
#' @noRd
.disk_raw_detour <- function(cell_area, n_valid) {
    ref_valid <- .reference_disk_raster(cell_area, n_valid)
    ref_area <- as.numeric(terra::global(ref_valid, "sum", na.rm = TRUE)[1, 1]) * cell_area
    ref_hull <- terra::convHull(.valid_polygons(ref_valid))
    ref_hull_perimeter <- .muffle_crs_warnings(terra::perim(ref_hull))
    (2 * sqrt(pi * ref_area)) / ref_hull_perimeter
}

#' `gm_reock_index()`'s reference: the raw `area/mbc_area` score of an
#' actual rasterized disk at the input's own cell size and area - a single
#' realization, unlike `.phase_averaged_ref_pp()` above (verified low
#' phase sensitivity, see file header).
#' @param cell_area scalar
#' @param n_valid target valid-cell count for the reference disk
#' @return raw reock score of the reference disk
#' @noRd
.disk_raw_reock <- function(cell_area, n_valid) {
    ref_valid <- .reference_disk_raster(cell_area, n_valid)
    ref_area <- as.numeric(terra::global(ref_valid, "sum", na.rm = TRUE)[1, 1]) * cell_area
    ref_hull <- terra::convHull(.valid_polygons(ref_valid))
    ref_pts <- unique(terra::crds(ref_hull))
    ref_mec <- .min_enclosing_circle(ref_pts)
    ref_area / (pi * ref_mec$r^2)
}

## -- the six scores ---------------------------------------------------------

#' Convex hull area ratio of a terra raster
#'
#' area / area(convex hull), in `(0, 1]`, `1` = the shape is convex (equals
#' its own convex hull). No `weighted` argument - see this file's own
#' header for why area/hull are structurally unweighted quantities.
#' @param rast a terra SpatRaster. The shape is derived directly from
#'   `rast`: a cell is part of the shape iff its own value is neither `NA`
#'   nor exactly `0` - both are holes, no separate mask argument.
#' @return `list(index, area, hull_area, hull, n_valid_cells)`. `hull` is
#'   the convex hull itself, a terra SpatVector, for plotting.
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_hull_ratio_index(r)$index
#' @export
gm_hull_ratio_index <- function(rast) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_hull_ratio_index")
    valid <- .valid_cells(rast)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, area = 0, hull_area = NA_real_, hull = NULL, n_valid_cells = 0L))
    }

    hull <- terra::convHull(.valid_polygons(valid))
    hull_area <- .muffle_crs_warnings(terra::expanse(hull))
    index <- if (hull_area > 0) area / hull_area else NA_real_
    list(index = index, area = area, hull_area = hull_area, hull = hull, n_valid_cells = as.integer(n_valid))
}

#' Polsby-Popper compactness score of a terra raster
#'
#' `raw / ref_index`, in `(0, 1]`, `1` = a circle, where `raw` is the
#' textbook `4*pi*area/perimeter^2` (perimeter via marching squares,
#' `terra::as.contour()` - see this file's own header for why, verified
#' empirically against the naive cell-edge alternative) and `ref_index` is
#' that SAME formula's own score on an actual rasterized disk at `rast`'s
#' own cell size and area.
#'
#' WHY NOT JUST `raw`: comparing it against the textbook constant `1`
#' (a claim about a perfect, infinitely-smooth circle) systematically
#' understates how compact a rasterized shape is. A rasterized DISK - the
#' one shape this score should never meaningfully penalize - scored `0.88`
#' on `raw` alone at a representative resolution, not because the shape is
#' imperfect, but because marching-squares perimeter overestimates a
#' raster boundary's true length, and squaring it in the denominator
#' roughly doubles that percentage error. Dividing by `ref_index` instead
#' cancels that measurement bias, since both numerator and reference go
#' through the identical perimeter-measurement pipeline: the same disk
#' scores `0.97-1.00` this way, across resolutions. This reference is
#' PHASE-AVERAGED (25 sub-cell offsets) rather than a single realization -
#' see this file's own header for why marching-squares perimeter needs
#' that where the other resolution-matched references in this file don't.
#' No `analytical_ref` opt-out is provided: the raw closed-form score is
#' not a more "correct" number to fall back to, it is measurably the wrong
#' one for a rasterized shape at any finite resolution.
#' No `weighted` argument - see this file's own header for why.
#' @inheritParams gm_hull_ratio_index
#' @return `list(index, area, perimeter, ref_index, n_valid_cells)`.
#'   `ref_index` is the reference disk's own raw score - `index` itself is
#'   `(4*pi*area/perimeter^2) / ref_index`, not the raw formula directly.
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_polsby_popper_index(r)$index
#' @export
gm_polsby_popper_index <- function(rast) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_polsby_popper_index")
    valid <- .valid_cells(rast)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, area = 0, perimeter = NA_real_, ref_index = NA_real_, n_valid_cells = 0L))
    }

    perimeter <- .contour_perimeter(valid)
    raw <- if (perimeter > 0) (4 * pi * area) / perimeter^2 else NA_real_
    ref_index <- if (is.na(raw)) NA_real_ else .phase_averaged_ref_pp(cell_area, n_valid)
    index <- if (is.na(raw)) NA_real_ else raw / ref_index
    list(index = index, area = area, perimeter = perimeter, ref_index = ref_index, n_valid_cells = as.integer(n_valid))
}

#' Width-length ratio of a terra raster's bounding box
#'
#' The shorter of the bounding box's x/y extents over the longer, in
#' `(0, 1]`, `1` = a square bounding box. Axis-aligned, not the minimum
#' bounding rectangle at any rotation - a diagonally-oriented elongated
#' shape can score deceptively high, the classic limitation of this score.
#' No `weighted` argument - see this file's own header for why.
#'
#' KNOWN LIMITATION (ported unchanged from shapeindices): blind to both
#' holes and multi-part dispersal, since only the bounding box's own
#' extent enters the ratio - a shape and the same shape with a large hole
#' punched through it score identically.
#' @inheritParams gm_hull_ratio_index
#' @return `list(index, length, width, n_valid_cells)`
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_width_length_ratio_index(r)$index
#' @export
gm_width_length_ratio_index <- function(rast) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_width_length_ratio_index")
    valid <- .valid_cells(rast)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, length = NA_real_, width = NA_real_, n_valid_cells = 0L))
    }

    e <- terra::ext(terra::trim(terra::ifel(valid, 1, NA)))
    dx <- e$xmax - e$xmin
    dy <- e$ymax - e$ymin
    width <- min(dx, dy)
    length <- max(dx, dy)
    index <- if (length > 0) width / length else NA_real_
    list(index = index, length = length, width = width, n_valid_cells = as.integer(n_valid))
}

#' Reock compactness score of a terra raster
#'
#' `raw / ref_index`, in `(0, 1]`, `1` = a circle, where `raw` is
#' `area/area(minimum bounding circle)` and `ref_index` is that SAME
#' formula's own score on an actual rasterized disk at `rast`'s own cell
#' size and area. The minimum bounding circle depends only on the shape's
#' own convex hull vertices, found via a deterministic (no RNG) Welzl-style
#' algorithm - exact, not an approximation.
#'
#' WHY NOT JUST `raw`: `area` is a pixel COUNT, which undershoots a
#' shape's true continuous area near its own boundary (a boundary cell
#' only counts as "in" when its centre does), while `mbc_area` comes from
#' EXACT hull vertices - two different measurement conventions for what
#' should be the same disk. A rasterized DISK scored `raw = 0.953` at a
#' representative resolution, not because the shape is imperfect, but from
#' that convention mismatch alone. Dividing by `ref_index` instead cancels
#' it, since both numerator and reference go through the identical
#' pixel-count-area / exact-MBC-area pipeline: the same disk scores
#' `0.996-1.000` this way, across resolutions (verified: low phase
#' sensitivity here, unlike `gm_polsby_popper_index()`'s own marching-
#' squares perimeter - see this file's own header). No `analytical_ref`
#' opt-out is provided: the raw ratio is not a more "correct" number to
#' fall back to, it is measurably the wrong one for a rasterized shape at
#' any finite resolution.
#' No `weighted` argument - see this file's own header for why.
#' @inheritParams gm_hull_ratio_index
#' @return `list(index, area, mbc_area, ref_index, mbc, n_valid_cells)`.
#'   `mbc` is the minimum bounding circle itself, a terra SpatVector, for
#'   plotting. `ref_index` is the reference disk's own raw score - `index`
#'   itself is `(area/mbc_area) / ref_index`, not the raw ratio directly.
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_reock_index(r)$index
#' @export
gm_reock_index <- function(rast) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_reock_index")
    valid <- .valid_cells(rast)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, area = 0, mbc_area = NA_real_, ref_index = NA_real_, mbc = NULL, n_valid_cells = 0L))
    }

    hull <- terra::convHull(.valid_polygons(valid))
    pts <- unique(terra::crds(hull))
    mec <- .min_enclosing_circle(pts)
    mbc_area <- pi * mec$r^2
    raw <- if (mbc_area > 0) area / mbc_area else NA_real_
    ref_index <- if (is.na(raw)) NA_real_ else .disk_raw_reock(cell_area, n_valid)
    index <- if (is.na(raw)) NA_real_ else raw / ref_index
    mbc_center <- terra::vect(matrix(mec$center, nrow = 1), type = "points", crs = terra::crs(rast))
    mbc <- terra::buffer(mbc_center, width = mec$r)
    list(index = index, area = area, mbc_area = mbc_area, ref_index = ref_index, mbc = mbc, n_valid_cells = as.integer(n_valid))
}

#' Detour compactness score of a terra raster
#'
#' `raw / ref_index`, in `(0, 1]`, `1` = a circle, where `raw` is the
#' ratio of the equal-area circle's ANALYTICAL perimeter to the perimeter
#' of the shape's own convex hull, and `ref_index` is that SAME formula's
#' own score on an actual rasterized disk at `rast`'s own cell size and
#' area. Angel, Parent & Civco (2010) introduce this to measure how hard a
#' shape is to circumnavigate as an obstacle. The hull's own perimeter
#' needs no marching-squares correction (unlike
#' `gm_polsby_popper_index()`'s raw shape perimeter) - `terra::convHull()`'s
#' output is already a proper straight-edged polygon, not a raster-traced
#' boundary.
#'
#' WHY NOT JUST `raw`: under centre-based rasterization, a boundary cell's
#' own far CORNER can sit outside the shape's true continuous boundary by
#' up to half a cell's diagonal, so the hull built from those corners has
#' vertices that genuinely protrude past the true circle - its perimeter
#' comes out LARGER than the analytical circle's own, even for an actual
#' disk. That scored `raw = 0.982` at a representative resolution.
#' Dividing by `ref_index` instead cancels it, since both numerator and
#' reference go through the identical hull-perimeter pipeline: the same
#' disk scores `0.997-1.000` this way, across resolutions (verified: low
#' phase sensitivity here, unlike `gm_polsby_popper_index()`'s own
#' marching-squares perimeter - see this file's own header). No
#' `analytical_ref` opt-out is provided: the raw ratio is not a more
#' "correct" number to fall back to, it is measurably the wrong one for a
#' rasterized shape at any finite resolution.
#' No `weighted` argument - see this file's own header for why.
#' @inheritParams gm_hull_ratio_index
#' @return `list(index, area, hull_perimeter, ref_index, hull, n_valid_cells)`.
#'   `hull` is the convex hull itself, a terra SpatVector, for plotting.
#'   `ref_index` is the reference disk's own raw score - `index` itself is
#'   `(2*sqrt(pi*area)/hull_perimeter) / ref_index`, not the raw ratio
#'   directly.
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_detour_index(r)$index
#' @export
gm_detour_index <- function(rast) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_detour_index")
    valid <- .valid_cells(rast)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, area = 0, hull_perimeter = NA_real_, ref_index = NA_real_, hull = NULL, n_valid_cells = 0L))
    }

    hull <- terra::convHull(.valid_polygons(valid))
    hull_perimeter <- .muffle_crs_warnings(terra::perim(hull))
    circle_perimeter <- 2 * sqrt(pi * area)
    raw <- if (hull_perimeter > 0) circle_perimeter / hull_perimeter else NA_real_
    ref_index <- if (is.na(raw)) NA_real_ else .disk_raw_detour(cell_area, n_valid)
    index <- if (is.na(raw)) NA_real_ else raw / ref_index
    list(index = index, area = area, hull_perimeter = hull_perimeter, ref_index = ref_index, hull = hull, n_valid_cells = as.integer(n_valid))
}

#' Exchange compactness score of a terra raster
#'
#' The share of the shape's own area that falls inside the equal-area
#' circle centred at its (unweighted) centroid, in `[0, 1]`, `1` = a
#' circle. Angel, Parent & Civco (2010) introduce this as a natural
#' metric for gerrymandering. Computed via pixel counting - each valid
#' cell's own centre is tested against the reference circle directly, the
#' raster-native analogue of the vector package's exact polygon
#' intersection. No `weighted` argument - see this file's own header for
#' why.
#'
#' KNOWN LIMITATION (ported unchanged from shapeindices): for a multi-part
#' shape, the reference circle is centred at the OVERALL centroid, which
#' for well-separated parts can sit in the empty space between them - once
#' the circle's radius is smaller than the distance to the nearest part,
#' the index is exactly 0, not just low.
#' @inheritParams gm_hull_ratio_index
#' @return `list(index, area, circle_area, circle, n_valid_cells)`.
#'   `circle` is the equal-area reference circle itself, a terra
#'   SpatVector, for plotting.
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_exchange_index(r)$index
#' @export
gm_exchange_index <- function(rast) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_exchange_index")
    valid <- .valid_cells(rast)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, area = 0, circle_area = NA_real_, circle = NULL, n_valid_cells = 0L))
    }

    cc <- .coord_rasters(valid)
    Gx <- as.numeric(terra::global(terra::ifel(valid, cc$x, NA), "mean", na.rm = TRUE)[1, 1])
    Gy <- as.numeric(terra::global(terra::ifel(valid, cc$y, NA), "mean", na.rm = TRUE)[1, 1])

    r_eq <- sqrt(area / pi)
    dist_to_g <- sqrt((cc$x - Gx)^2 + (cc$y - Gy)^2)
    inside <- valid & (dist_to_g <= r_eq)
    n_inside <- as.numeric(terra::global(inside, "sum", na.rm = TRUE)[1, 1])
    inter_area <- n_inside * cell_area
    index <- inter_area / area

    circle_center <- terra::vect(matrix(c(Gx, Gy), nrow = 1), type = "points", crs = terra::crs(rast))
    circle <- terra::buffer(circle_center, width = r_eq)
    list(index = index, area = area, circle_area = pi * r_eq^2, circle = circle, n_valid_cells = as.integer(n_valid))
}
