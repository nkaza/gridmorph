## =========================================================================
## Depth index: mean distance to the boundary vs. an equal-area circle,
## ported from shapeindices::depth_index(). index = mean_depth / ref_depth,
## in (0, 1] - actual/ref, not ref/actual, because the disk MAXIMISES this
## quantity (Angel, Parent & Civco 2010's Depth Proposition), the opposite
## direction from every other index in this package.
##
## gm_depth_index() is built directly on this file's own distance-transform
## primitive rather than the morphology module's erode()/dilate() - the
## raw distance transform is the field itself, not a derived erosion/
## dilation result. Uses terra::distance() (chunk-aware on disk-backed
## rasters, no mmand dependency needed) - cross-checked against
## mmand::distanceTransform(signed = TRUE) on a hand-checkable test case
## before trusting it: max absolute difference exactly 0. Same "prefer
## terra-native, chunk-safe operations over mmand's in-memory-array
## requirement" decision as R/morphology.R's own erode()/dilate() switch -
## see that file's header for the full reasoning and the box/disc-kernel
## speed comparison that didn't change the conclusion.
##
## RASTER-EDGE PADDING: `.depth_field()` (below) pads before calling
## terra::distance() - a real, verified bug otherwise, the distance()
## analogue of R/morphology.R's own focal() fillvalue = 0 fix (a
## different terra function, same underlying problem: a shape reaching
## the raster's true edge has no boundary there for distance() to measure
## against, understating its own depth - verified directly against a
## value outside this index's own documented (0, 1] range). See
## .depth_field()'s own comment for the full reproduction.
##
## REFERENCE IS RESOLUTION-MATCHED, NOT CLOSED-FORM - the one place this
## index's reference construction genuinely differs from every other index
## in this package (moment/span/radial-concentration all use a pure
## closed-form or exact-quadrature reference). Reason: `mean_depth` itself
## comes from `.depth_field()`'s `terra::distance()` cell-center Euclidean
## distance transform, which carries a real, systematic, resolution-
## convergent bias relative to true continuous distance-to-boundary
## (confirmed by decomposition: distance-to-nearest-invalid-CELL-CENTRE
## overstates distance-to-the-cell's-own-EDGE by `~0.6 * cell_size`
## regardless of depth, while the rasterized boundary itself sits inside
## the true continuous boundary by `~0.2-0.4 * cell_size` more at greater
## depth - both effects genuinely O(cell_size), confirmed by holding
## physical radius fixed and varying cell_size alone). Comparing that
## biased numerator against an unbiased closed-form reference does not
## cancel the bias - it leaked straight into the index (a rasterized disk
## scored `1.024`, not `1`, at the resolution `vignette("a-basic-usage")`
## uses). A CONSTANT correction doesn't work either (tried and rejected:
## subtracting `0.5 * cell_size` from `mean_depth` moves the same disk to
## `0.97`, since the two error terms partly cancel each other already and
## a flat subtraction only removes one of them).
##
## Fix: build an ACTUAL disk raster at the input's own cell size and area,
## arrange its cells as concentric density-sorted rings (innermost =
## highest density, the same optimal arrangement the closed-form reference
## assumes), and run it through the EXACT SAME `.depth_field()` call the
## input itself uses (`.resolution_matched_ref_depth()`, below). Both
## `mean_depth` and `ref_depth` then carry the same systematic bias, so it
## cancels in their ratio instead of leaking into `index` - verified
## directly: the same disk that scored `1.024` against the closed-form
## reference scores `1.001` against this one, and the improvement holds
## across a 4x resolution range (`1.011 -> 1.001`, `1.004 -> 1.0001`) and
## is insensitive to the reference disk's own sub-cell centring (`+-0.001`
## across widely different phases) - not a coincidence specific to one
## test case. `gm_radial_concentration_index()` does NOT need this: its
## own `D1_ref` (`.annulus_reference_D1_raster()`, R/utils.R) is compared
## against a point-sampled `D1` that carries no analogous distance-
## transform bias, so the plain closed form is already correctly unbiased
## there and is left untouched.
##
## Every valid cell has the SAME area (unlike triangles), which simplifies
## the ring construction the reference raster's cells get sorted into:
## sorting N cells by density descending, the k-th ring's outer radius is
## exactly `sqrt(k * cell_area / pi)`. ADAPTIVE (see R/utils.R's
## `.adaptive_density_bins()`): exact per-cell sort when the raster is
## small enough that sorting is cheap, falling back to a memory-bounded
## K-bin histogram approximation only above that threshold - never pulls a
## full per-cell vector into memory for a genuinely large raster.
##
## ALWAYS the "weighted" computation, `weighted = FALSE` included: passing
## a constant-1 density (see `.mass_raster()`, R/utils.R) collapses the
## ring assignment to a single ring spanning the whole reference disk, so
## `ref_depth` reduces to that disk's own plain (unweighted) mean depth -
## no separate closed-form code path to maintain alongside this one. This
## no longer reduces to the exact closed form `R/3` algebraically the way
## the previous closed-form reference did (it's a raster computation now,
## not a formula) - it converges to `R/3` at fine resolution, same as
## `mean_depth` itself does, which is the entire point.

#' Reference mean depth, computed on an ACTUAL rasterized disk at the same
#' cell size and area as the input, run through the same `.depth_field()`
#' distance transform - not a closed-form integral. See file header for
#' why: `.depth_field()` carries a real, systematic, resolution-convergent
#' bias, and only a reference built from the SAME biased pipeline cancels
#' it in `mean_depth / ref_depth`.
#'
#' The disk's valid cells are assigned density by cumulative FRACTION of
#' area (innermost = `density_sorted[1]`), not absolute rank - rasterizing
#' a disk at a given radius essentially never lands on exactly `n_total`
#' cells, so a fraction-based lookup reproportions ring widths slightly
#' instead of breaking the ring assignment outright when the reference
#' raster's own cell count differs from `n_total`.
#' @param cell_area scalar
#' @param density_sorted numeric vector, density per ring, sorted
#'   descending
#' @param count numeric vector, number of cells in each ring (default: one
#'   cell per ring)
#' @return ref_depth
#' @noRd
.resolution_matched_ref_depth <- function(cell_area, density_sorted, count = rep(1, length(density_sorted))) {
    n_total <- sum(count)
    ref_valid <- .reference_disk_raster(cell_area, n_total)

    depth <- .depth_field(ref_valid)
    valid_v <- as.vector(terra::values(ref_valid))
    depth_v <- as.vector(terra::values(depth))[valid_v]
    n_ref <- length(depth_v)

    cum_frac <- cumsum(count) / n_total
    frac <- (rank(-depth_v, ties.method = "first") - 0.5) / n_ref
    ring_idx <- pmin(findInterval(frac, cum_frac) + 1L, length(density_sorted))
    w <- density_sorted[ring_idx]

    sum(w * depth_v) / sum(w)
}

#' Euclidean distance transform of `valid` (TRUE = inside), in the physical
#' units of `rast`'s own CRS, via `terra::distance()` on an inverted mask -
#' `terra::distance()` computes distance FROM `NA` cells TO the nearest
#' non-`NA` cell, so invert (valid -> `NA`, invalid -> a real marker value)
#' to get distance FROM inside cells TO the boundary. Chunk-aware on
#' disk-backed rasters the same way `focal()` is (see `R/morphology.R`'s own
#' header for why that matters more than raw per-call speed for this
#' package) - no `mmand` dependency needed here either. Verified against
#' `mmand::distanceTransform(signed = TRUE)` on a hand-checked test case
#' before switching: max absolute difference exactly 0.
#'
#' PADS `inverted` by one cell (`terra::extend(..., fill = 1)`) before
#' calling `distance()`, then crops back - a real, verified bug otherwise,
#' the `terra::distance()` analogue of `R/morphology.R`'s own `focal()`
#' `fillvalue` fix (a different terra function, same underlying problem:
#' `distance()` has no `fillvalue`-equivalent argument at all, so a shape
#' reaching `valid`'s own raster edge has no boundary to measure distance
#' to there, understating its own depth). Verified directly: an edge-
#' touching square gave `depth_index = 1.238` - a value outside this
#' index's own documented `(0, 1]` range, since a square must always score
#' below the disk-optimal `1` - while the padded fix gives EXACTLY the
#' same `mean_depth` as the identical square positioned away from the
#' edge (translation invariance, which any correct index must have). Also
#' fixes a related degenerate case: a shape that fills its ENTIRE raster
#' (no invalid cell anywhere for `distance()` to measure to) previously
#' returned `NaN`, now correctly measures against the padding instead.
#' @param valid logical terra SpatRaster (TRUE = inside the shape)
#' @return terra SpatRaster, same grid as `valid`, values = depth (>= 0) at
#'   TRUE cells, `NA` elsewhere
#' @noRd
.depth_field <- function(valid) {
    inverted <- terra::ifel(valid, NA, 1)
    padded <- terra::extend(inverted, 1, fill = 1)
    # .muffle_crs_warnings(): terra::distance() emits its own "unknown crs"
    # warning on a missing-CRS raster, duplicating the clearer one
    # .check_planar_crs() already gave the caller at the top of
    # gm_depth_index() - the same class of redundant-warning gap
    # .contour_perimeter() (classical-metrics.R) already guards against
    # for terra::perim().
    d <- terra::crop(.muffle_crs_warnings(terra::distance(padded)), valid)
    terra::ifel(valid, d, NA)
}

#' Depth index for a terra raster
#'
#' Mean distance from inside points to the shape's own boundary, relative
#' to an equal-area circle. index in `(0, 1]`, `= 1` iff the shape is
#' (almost everywhere) a disk - see the file header for the
#' Brunn-Minkowski-based proof this ports from
#' `shapeindices::depth_index()`. Handles holes and multi-part shapes with
#' no special-casing at all: they're just cells with different values in
#' one grid.
#' @param rast a terra SpatRaster. The shape is derived directly from
#'   `rast`: a cell is part of the shape iff its own value is neither `NA`
#'   nor exactly `0` - both are holes, no separate mask argument.
#' @param weighted logical, default `TRUE`. When `TRUE`, `rast`'s own cell
#'   values are used as the density/mass throughout. When `FALSE`, every
#'   valid cell is treated as equally massed regardless of `rast`'s actual
#'   values - exactly reproducing the plain (unweighted) index even if
#'   `rast` is itself a continuous raster. If `rast` is a CATEGORICAL
#'   raster (`terra::is.factor(rast)` is `TRUE` - e.g. a land-use
#'   classification with an attached levels table) and `weighted = TRUE`,
#'   the category CODE NUMBERS themselves get used as mass, which is
#'   rarely the intended comparison (a class coded `3` outweighs one coded
#'   `1` for no meaningful reason) - a warning is issued but the
#'   computation still proceeds on the literal codes, since `weighted =
#'   TRUE` was explicitly requested; pass `weighted = FALSE` for
#'   categorical rasters unless the codes really are meant as magnitudes.
#'   This check never fires on non-categorical (including plain binary
#'   0/1) rasters, so `weighted = TRUE` is always safe there.
#' @param n_bins integer, the exact/binned threshold and bin count for the
#'   weighted reference's concentric-rings construction (only matters when
#'   `weighted = TRUE` and the density genuinely varies - `weighted = FALSE`
#'   always collapses to a single ring regardless of `n_bins`) - exact (no
#'   accuracy cost) when the number of valid cells is `<= n_bins`, a K-bin
#'   histogram approximation above that. Default `1000` measured at ~0.1%
#'   relative error on a genuinely non-uniform weight pattern; see
#'   `.adaptive_density_bins()` in `R/utils.R`.
#' @return list(index, mean_depth, ref_depth, area, n_valid_cells). Computing
#'   `ref_depth` costs one more `terra::distance()` call, on a freshly
#'   rasterized disk of the same cell size and area as `rast` (not on
#'   `rast` itself) - see file header for why a closed-form reference isn't
#'   used here the way every other index in this package uses one.
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_depth_index(r)$index
#' @export
gm_depth_index <- function(rast, weighted = TRUE, n_bins = 1000) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_planar_crs(rast, "gm_depth_index")
    valid <- .valid_cells(rast)
    w <- .mass_raster(rast, valid, weighted)

    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    cell_area <- prod(terra::res(rast))
    area <- n_valid * cell_area

    if (n_valid == 0) {
        warning("No valid cells; index is not defined.")
        return(list(index = NA_real_, mean_depth = NA_real_, ref_depth = NA_real_,
                    area = 0, n_valid_cells = 0L))
    }

    depth <- .depth_field(valid)
    valid_v <- as.vector(terra::values(valid))
    depth_v <- as.vector(terra::values(depth))[valid_v]

    w_v <- as.vector(terra::values(w))[valid_v]
    mean_depth <- sum(w_v * depth_v) / sum(w_v)
    bins <- .adaptive_density_bins(valid, w, n_bins)
    ref_depth <- .resolution_matched_ref_depth(cell_area, bins$density, bins$count)

    list(index = mean_depth / ref_depth, mean_depth = mean_depth, ref_depth = ref_depth,
         area = area, n_valid_cells = as.integer(n_valid))
}
