## =========================================================================
## Shared helpers: valid-cell derivation and temp-file hygiene. Every
## exported index function routes through these rather than re-implementing
## its own version of the same checks.
## =========================================================================
##
## SINGLE-RASTER INPUT: every index function takes ONE `rast` argument, not
## a separate binary mask plus an optional weight overlay. The shape itself
## is DERIVED from `rast`: a cell is part of the shape iff its value is
## neither `NA` nor exactly `0` - `NA` and `0` are both holes, the same
## precedent `shapeindices` already documents for its own (separate)
## `weight` argument there.
##
## `weighted = TRUE/FALSE` (every exported index function's own second
## argument) then controls what MASS each valid cell carries, not which
## cells are valid - that part never changes: `weighted = TRUE` (default)
## uses `rast`'s own cell values as the density/mass throughout;
## `weighted = FALSE` ignores their magnitude and treats every valid cell
## as equally massed, exactly reproducing the OLD unweighted calling
## convention (a separate binary mask, no weight argument at all) - not
## approximately: every index's own weighted computation is provably exact
## for a spatially-constant density (see `.bin_index_raster()`'s own
## "degenerate" mode, and the file-level comment in each index file working
## through the algebra), so there's no separate closed-form code path left
## to maintain - `weighted = FALSE` reuses the exact same machinery with a
## constant-1 density raster in place of `rast`'s own values.
##
## This also means no cross-raster grid alignment is needed at all - there
## is only ever one raster now, nothing to register against anything else.
## (An earlier design had a separate mask + weight, with automatic
## resampling/CRS-checking between them - removed entirely along with the
## second raster it existed to reconcile.)
##
## PLANAR CRS REQUIRED (FOR NOW): NOT because area, distance, or
## perimeter can't be computed correctly on geographic (longitude/
## latitude) data - terra's OWN primitives (`terra::distance()`,
## `terra::perim()`, `terra::expanse()`) compute all three correctly on a
## genuine geographic CRS, via geodesic (ellipsoid-aware) math, not naive
## degree arithmetic - this is terra's own documented, deliberate
## behaviour (see `?terra::expanse`: "the best way to compute area is to
## use the longitude/latitude crs if that is what the data come in...
## contrary to (erroneous) popular belief [that you should reproject to a
## planar CRS first]").
##
## THE REAL PROBLEM: line-crossing tests and point-to-point distance
## accumulation done in RAW (x, y) coordinate space, never routed through
## any of terra's own geodesic-aware primitives at all.
## `gm_convexity_index()` tests whether a straight LINE between two
## sampled points stays inside the shape by linearly interpolating their
## raw x/y coordinates (`.frac_outside_line()`) - a straight line in
## (longitude, latitude) space is not a geodesic (great-circle) path, so
## this test means something different, and wrong, on a sphere.
## `gm_span_index()`/`gm_radial_concentration_index()`'s own pairwise
## distances (`sqrt(dx^2 + dy^2)` on raw coordinates), moment-indices.R's
## Ixx/Iyy/Ixy computation (built from raw coordinate DIFFERENCES about
## the centroid), and classical-metrics.R's minimum-enclosing-circle fit
## (`.circle_from_2()`/`.circle_from_3()`, pure planar circle geometry)
## are all the same kind of point-based Euclidean math, assumed
## throughout, with no geodesic awareness at all. Even
## `gm_depth_index()`'s OWN area computation (`cell_area =
## prod(terra::res(rast))`, nowhere near `terra::distance()` but still
## combined with its own geodesic-metres output) is naive degree
## arithmetic, not `terra::expanse()`'s correct geodesic area - a real,
## verified bug (index = 107695.9 for a value documented to be in
## `(0, 1]`), but the actual fix for THAT specific case would be "use
## terra's own area function", not "area is impossible in degrees".
##
## `.check_planar_crs()` (below) still hard-errors on a genuine geographic
## CRS, for every exported index alike, because the point-based Euclidean
## assumptions above run through essentially every index in this package
## in one form or another, and untangling which ones could be salvaged
## with a geodesic-aware rewrite (some plausibly could - `gm_depth_index()`
## itself might need only its own area formula fixed; convexity's line-
## crossing test and the MEC circle-fitting math are the genuinely hard
## cases) is real, not-yet-attempted future work. Reproject to a planar
## CRS first for now; that's not a workaround for a limitation terra
## itself also has, it's THIS package's own, not-yet-closed gap.
##
## MISSING CRS (`terra::crs(rast) == ""`) is different: `terra::is.lonlat()`
## returns `NA` for it (not `TRUE`), and verified directly that every
## terra primitive this package uses (`distance()`, `perim()`, `expanse()`)
## already falls back to ordinary Euclidean math in that case - correct,
## not silently wrong, for the common case of a genuinely abstract/
## synthetic grid with no real-world geographic meaning (exactly what
## `crs = "local"` is for, used throughout this package's own tests and
## examples). Warned about, not hard-errored: unlike a real geographic
## CRS, there is no computation here that's actually WRONG, only
## unverifiable - only a `NA`/missing CRS could ALSO mean "this raster's
## own CRS was lost somewhere upstream and it really is geographic data",
## which the caller is in a better position to know than this package is.

#' Evaluate `expr`, muffling any warning whose own message mentions
#' "crs" (case-insensitively) - used to swallow terra's own scattered,
#' inconsistent "unknown crs" chatter (from `is.lonlat()`, `expanse()`,
#' `perim()` - each warns separately, in its own wording, if called on a
#' missing-CRS raster) once `.check_planar_crs()` has already told the
#' caller about the SAME missing-CRS situation clearly, at the top of the
#' function. Message-matched, not a blanket `suppressWarnings()`: any
#' OTHER, unrelated warning from the wrapped expression still comes
#' through normally.
#' @param expr an expression to evaluate
#' @return the value of `expr`
#' @noRd
.muffle_crs_warnings <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
        if (grepl("crs", conditionMessage(w), ignore.case = TRUE)) {
            invokeRestart("muffleWarning")
        }
    })
}

#' Hard-error on a genuine geographic (longitude/latitude) CRS. NOT because
#' area/distance/perimeter can't be computed correctly there - terra's own
#' `distance()`/`perim()`/`expanse()` do that correctly via geodesic math,
#' confirmed against terra's own documentation - but because this
#' package's own point-based Euclidean geometry (line-crossing tests,
#' pairwise distance accumulation, moment-tensor math, circle-fitting)
#' assumes raw x/y values are planar Cartesian throughout, and is not
#' geodesic-aware (see this file's own header for exactly which
#' computations and the two distinct, verified failure modes this closes
#' off). Warns, rather than erroring, if `rast` has no CRS at all -
#' verified that every terra primitive this package relies on already
#' computes correctly (plain Euclidean) in that case, so there's nothing
#' actually WRONG to block, only something worth flagging in case the
#' missing CRS wasn't intentional.
#' @param rast terra SpatRaster, the exported function's own primary input
#' @param fn_name character, the calling function's own name, for the
#'   error/warning message
#' @return invisible `NULL`
#' @noRd
.check_planar_crs <- function(rast, fn_name) {
    # terra::is.lonlat()'s own `warn` parameter does NOT suppress this
    # particular "unknown crs" warning - verified directly by reading its
    # source: that warning fires unconditionally in the missing-CRS
    # branch, only the DIFFERENT "assuming lon/lat crs" warning (in its
    # unrelated `perhaps = TRUE` heuristic-guessing mode, never used here)
    # is actually gated by `warn`. `.muffle_crs_warnings()` below gets the
    # single, clearer gridmorph-level warning through without terra's own
    # redundant duplicate right alongside it.
    lonlat <- .muffle_crs_warnings(terra::is.lonlat(rast, warn = FALSE))
    if (isTRUE(lonlat)) {
        stop(fn_name, "(): `rast` has a geographic (longitude/latitude) CRS (", terra::crs(rast, describe = TRUE)$name,
             "). This is not because area, distance, or perimeter can't be computed correctly on geographic ",
             "data - terra's own distance()/perim()/expanse() do that correctly, via geodesic math. The actual ",
             "problem is that this package's own point-based geometry (line-crossing tests for convexity, ",
             "pairwise distances for span/radial_concentration, moment-tensor calculations, minimum-enclosing- ",
             "circle fitting) works directly in raw x/y coordinate space throughout, which assumes planar ",
             "(Cartesian) coordinates and is not yet geodesic-aware. Reproject to a suitable projected CRS ",
             "first (e.g. via terra::project()) before calling this function - not silently corrected here.")
    }
    if (is.na(lonlat)) {
        warning(fn_name, "(): `rast` has no CRS set. Every computation here assumes `rast`'s own x/y ",
                "coordinate units already ARE physical planar distance (as they would be for a purely ",
                "synthetic/abstract grid) - proceeding on that assumption. If this raster represents ",
                "real-world geographic data, set its CRS explicitly (e.g. terra::crs(rast) <- \"EPSG:...\") ",
                "so a geographic CRS can be caught and rejected here instead of silently assumed planar.",
                call. = FALSE)
    }
    invisible(NULL)
}

#' The single boolean "valid cell" indicator every index computes over -
#' derived directly from `x`: part of the shape iff its value is neither
#' `NA` nor exactly `0`. One source of truth rather than ad hoc
#' `is.na()`/`== 0` checks scattered through each index.
#'
#' On a categorical (`terra::is.factor()`) `x`, this treats a cell coded
#' `0` as a hole even if `0` is one of the raster's own real class levels
#' (e.g. "0 = water") - deliberately unchanged for categorical input,
#' since every index's own notion of "shape" is built on this exact rule.
#' Different from R/morphology.R's own categorical convention (`NA` is the
#' ONLY background there; every level, including `0`, persists) - not an
#' inconsistency to reconcile, since the two modules answer different
#' questions (morphology.R reshapes a raster's own labels; this defines
#' what "the shape" even is for an index to score) and predate the
#' categorical-morphology design by this package's entire index-function
#' half, where "0 = hole" has been the load-bearing shape definition from
#' the start.
#' @param x terra SpatRaster (an exported function's own `rast` argument)
#' @return logical terra SpatRaster
#' @noRd
.valid_cells <- function(x) {
    if (!inherits(x, "SpatRaster")) stop("`rast` must be a terra SpatRaster.")
    !is.na(x) & (x != 0)
}

#' Build an actual rasterized disk at a given cell size and (approximately)
#' target valid-cell count, centred on its own raster - shared by every
#' index whose own reference needs computing through the SAME raster
#' measurement pipeline as the input, rather than a pure closed-form
#' formula (see `R/depth-index.R`'s file header, and each classical-metric
#' function's own roxygen in `R/classical-metrics.R`, for why this matters
#' more than raw accuracy loss would suggest). Never lands on exactly
#' `n_valid` cells - rasterizing a circle at a given radius essentially
#' never does - but the resulting mismatch was verified to leave the
#' reference computation stable (see `R/depth-index.R`'s own header for the
#' sub-cell-phase sensitivity check this was confirmed against).
#' @param cell_area scalar, one cell's area
#' @param n_valid target valid-cell count (approximate)
#' @param phase length-2 numeric, `c(x, y)` offset of the disk's own centre
#'   from the raster's centre, in CELL-SIZE units. Default `c(0, 0)`
#'   (the raster-centred disk every caller except
#'   `gm_polsby_popper_index()`'s own phase-averaged reference uses - see
#'   `R/classical-metrics.R`'s file header for why that one needs several
#'   different phases averaged together and this one doesn't)
#' @return logical terra SpatRaster (TRUE = inside the disk), the same
#'   convention `.valid_cells()` returns
#' @noRd
.reference_disk_raster <- function(cell_area, n_valid, phase = c(0, 0)) {
    cell_size <- sqrt(cell_area)
    R <- sqrt(n_valid * cell_area / pi)
    pad <- ceiling(R / cell_size) + 2L
    n <- 2L * pad + 1L
    ext <- n * cell_size
    ref <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = ext, ymin = 0, ymax = ext, crs = "local")
    xy <- terra::xyFromCell(ref, seq_len(terra::ncell(ref)))
    cx <- ext / 2 + phase[1] * cell_size
    cy <- ext / 2 + phase[2] * cell_size
    d <- sqrt((xy[, 1] - cx)^2 + (xy[, 2] - cy)^2)
    terra::values(ref) <- ifelse(d <= R, 1, NA)
    !is.na(ref)
}

#' Warn (not error) when `weighted = TRUE` is requested on a categorical
#' raster - see `.mass_raster()`'s own comments for why. Extracted into
#' its own function (rather than left inline in `.mass_raster()`) so
#' `gm_shape_indices()` can also call it exactly ONCE up front and muffle
#' the redundant repeats from each of the seven weighted-capable index
#' functions it calls internally - the same "warn once, not seven times"
#' treatment `.check_planar_crs()` already gets there for a missing CRS.
#' Message text deliberately unchanged by that muffling:
#' `.muffle_categorical_warnings()` (below) matches on it directly.
#' @param rast terra SpatRaster
#' @return invisible `NULL`
#' @noRd
.warn_categorical_weighted <- function(rast) {
    warning("`rast` is a categorical raster (terra::is.factor(rast) is TRUE) and weighted = TRUE: ",
            "its category CODE NUMBERS will be used as mass/density directly, not just its footprint. ",
            "This is rarely the intended comparison for categorical data - pass weighted = FALSE to use ",
            "only where rast is nonzero/non-NA, ignoring the category codes themselves. Proceeding with ",
            "the codes as mass, since weighted = TRUE was explicitly requested.", call. = FALSE)
    invisible(NULL)
}

#' Evaluate `expr`, muffling any warning whose own message mentions
#' "categorical" (case-insensitively) - the `.warn_categorical_weighted()`
#' analogue of `.muffle_crs_warnings()`, for the same reason: swallow the
#' redundant repeats once `gm_shape_indices()` has already told the
#' caller about the SAME categorical-raster situation clearly, up front.
#' @param expr an expression to evaluate
#' @return the value of `expr`
#' @noRd
.muffle_categorical_warnings <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
        if (grepl("categorical", conditionMessage(w), ignore.case = TRUE)) {
            invokeRestart("muffleWarning")
        }
    })
}

#' The per-cell mass/density every index computes over, given the
#' `weighted` toggle every exported index function exposes: `rast`'s own
#' values when `weighted = TRUE`, or a constant 1 over every valid cell
#' when `FALSE` - see this file's own header for why the latter is an
#' EXACT, not approximate, reproduction of the old unweighted convention.
#'
#' `weighted = TRUE` on a raster `terra` itself considers categorical
#' (`terra::is.factor(rast)` - the raster carries an explicit `levels()`
#' table, e.g. a land-use classification read from a categorical GeoTIFF)
#' warns rather than erroring (via `.warn_categorical_weighted()`): the
#' category CODE NUMBERS get used as mass directly (verified: a two-class
#' raster coded 1 vs 3 gave a visibly different centroid/index than the
#' same shape with `weighted = FALSE`, 0.77 vs 0.96 on a real test case),
#' which is essentially never the intended comparison - but it's not
#' incoherent either (the numbers ARE real numbers, just not meant as a
#' magnitude), so the computation proceeds on the caller's explicit
#' `weighted = TRUE` request rather than silently overriding it.
#' `is.factor()` deliberately used over a value-distribution heuristic
#' (e.g. "few distinct integers"): it relies on the raster's own declared
#' type rather than guessing, so it never false-positives on legitimately
#' discrete weighted data (e.g. population counts 1, 2, 3, ...) that just
#' happens to have few distinct values.
#' @param rast terra SpatRaster, the exported function's own primary input
#' @param valid logical terra SpatRaster, `.valid_cells(rast)`
#' @param weighted logical
#' @return terra SpatRaster
#' @noRd
.mass_raster <- function(rast, valid, weighted) {
    if (!weighted) return(terra::ifel(valid, 1, NA))
    if (terra::is.factor(rast)) .warn_categorical_weighted(rast)
    rast
}

#' Adaptive exact/binned (density, count) pairs for a concentric-rings
#' weighted reference construction, shared by `gm_depth_index()`'s D1
#' reference and `gm_moment_of_inertia_index()`'s J_ref - both sort valid
#' cells by density descending and build rings from cumulative area,
#' differing only in what they do with the rings afterward (see each
#' file's own header).
#'
#' EXACT below `n_bins` valid cells (one ring per cell, `count = 1` each) -
#' zero accuracy cost for typical rasters, identical to always doing this
#' exactly. Above `n_bins`, bins `density_raster`'s own values into
#' `n_bins` levels via `terra::app()` + `terra::freq()` (chunk-safe, never
#' pulls the full per-cell vector into memory) and treats each populated
#' bin as one ring - verified empirically on a radial-gradient-plus-noise
#' test weight pattern: relative error vs. the exact computation is ~12%
#' at 10 bins, ~0.4% at 300, ~0.1% at 1000. Default `n_bins = 1000` trades
#' a real but small (~0.1%) accuracy cost, only above the threshold, for a
#' hard memory bound - not a hypothetical trade, a measured one.
#'
#' `terra::classify()` + `freq()` was tried and rejected for the binning
#' step: `freq()`'s `$value` column on a classified raster is a parsed
#' string label (`"[0 - 2]"`), and - the real trap - `freq()` silently
#' OMITS empty bins, so a naive positional bin-index mapping silently
#' misattributes counts whenever any bin happens to be empty (verified
#' directly: a 5-bin classification with a deliberately-empty middle bin
#' returned only 4 rows). Fixed by computing plain integer bin indices
#' directly via `terra::app()` instead - `freq()` on THAT raster gives a
#' numeric `$value` column trustable as the bin index itself, regardless
#' of which bins are empty.
#' Assign each valid cell to one of at most `n_bins` equal-width density
#' bins - the shared low-level construction behind `.adaptive_density_bins()`'s
#' (density, count) summary AND `.sample_valid_points()`'s per-bin
#' sampling (the latter needs the actual bin-index RASTER, not just the
#' aggregate summary, to draw real cell coordinates from a specific bin -
#' see that function's own header). Three mutually exclusive modes:
#'  - `"exact"`: `n_valid <= n_bins` - cheap to treat every cell as its
#'    own bin, no approximation at all
#'  - `"degenerate"`: density has zero range (e.g. a uniform weight) -
#'    every valid cell already IS one bin
#'  - `"binned"`: the general case - `r_idx` is a real bin-index raster
#' @param valid logical terra SpatRaster (TRUE = inside the shape)
#' @param density_raster terra SpatRaster, need not already be masked to
#'   `valid` - masked internally
#' @param n_bins integer, the exact/binned threshold and bin count
#' @return a mode-tagged list (see above)
#' @noRd
.bin_index_raster <- function(valid, density_raster, n_bins) {
    d <- terra::ifel(valid, density_raster, NA)
    n_valid <- as.numeric(terra::global(valid, "sum", na.rm = TRUE)[1, 1])
    if (n_valid <= n_bins) return(list(mode = "exact", d = d, n_valid = n_valid))

    rng <- as.numeric(terra::global(d, "range", na.rm = TRUE)[1, ])
    if (!all(is.finite(rng)) || rng[2] <= rng[1]) {
        return(list(mode = "degenerate", value = if (is.finite(rng[1])) rng[1] else 0, n_valid = n_valid))
    }

    bin_width <- (rng[2] - rng[1]) / n_bins
    r_idx <- terra::app(d, fun = function(x) pmin(floor((x - rng[1]) / bin_width), n_bins - 1))
    list(mode = "binned", r_idx = r_idx, rng = rng, bin_width = bin_width, n_valid = n_valid)
}

#' @param valid logical terra SpatRaster (TRUE = inside the shape)
#' @param density_raster terra SpatRaster, the SAME per-cell density
#'   quantity the caller's own moment/reference formula uses (e.g.
#'   `weight` itself for `gm_depth_index()`, or the normalised `rho` for
#'   `gm_moment_of_inertia_index()` - never recomputed/renormalised here).
#'   Need not already be masked to `valid` - masked internally.
#' @param n_bins integer, the exact/binned threshold and bin count
#' @return `list(density, count)`, both numeric vectors sorted
#'   density-descending, `sum(count) == n_valid`
#' @noRd
.adaptive_density_bins <- function(valid, density_raster, n_bins = 1000) {
    b <- .bin_index_raster(valid, density_raster, n_bins)

    if (b$mode == "exact") {
        d_v <- as.vector(terra::values(b$d))
        d_v <- d_v[!is.na(d_v)]
        return(list(density = sort(d_v, decreasing = TRUE), count = rep(1, length(d_v))))
    }
    if (b$mode == "degenerate") {
        return(list(density = b$value, count = b$n_valid))
    }

    f <- terra::freq(b$r_idx)
    bin_mid <- b$rng[1] + (f$value + 0.5) * b$bin_width
    ord <- order(bin_mid, decreasing = TRUE)
    list(density = bin_mid[ord], count = f$count[ord])
}

#' Best-effort estimate of currently-available system memory, in MB - used
#' only to size a Monte Carlo sample ceiling, never to guarantee anything,
#' so a wrong answer just makes the ceiling more/less conservative rather
#' than incorrect. Reads `/proc/meminfo`'s `MemAvailable` on Linux,
#' `vm_stat`'s free pages on macOS; falls back to a fixed conservative
#' constant everywhere else (including if either parse fails) - safe
#' because UNDER-estimating only shrinks the ceiling (more caution, never
#' unsafe), while OVER-estimating could reintroduce whatever memory
#' problem the caller is trying to avoid. Ported directly from
#' `shapeindices::.available_memory_mb()` - same OS-portability concern,
#' no `sf`/triangulation dependency to strip out.
#' @return numeric, estimated available memory in MB
#' @noRd
.available_memory_mb <- function() {
    fallback_mb <- 2048
    tryCatch({
        sys <- Sys.info()[["sysname"]]
        if (identical(sys, "Linux") && file.exists("/proc/meminfo")) {
            meminfo <- readLines("/proc/meminfo", n = 5L, warn = FALSE)
            avail <- grep("^MemAvailable:", meminfo, value = TRUE)
            kb <- if (length(avail) == 1) {
                as.numeric(regmatches(avail, regexpr("[0-9]+", avail)))
            } else NA_real_
            if (is.finite(kb) && kb > 0) return(kb / 1024)
            fallback_mb
        } else if (identical(sys, "Darwin")) {
            vmstat <- system("vm_stat", intern = TRUE, ignore.stderr = TRUE)
            page_line <- grep("page size of", vmstat, value = TRUE)
            page_bytes <- if (length(page_line) == 1) {
                as.numeric(regmatches(page_line, regexpr("[0-9]+", page_line)))
            } else NA_real_
            if (!is.finite(page_bytes) || page_bytes <= 0) page_bytes <- 4096
            free_line <- grep("^Pages free:", vmstat, value = TRUE)
            pages <- if (length(free_line) == 1) {
                as.numeric(regmatches(free_line, regexpr("[0-9]+", free_line)))
            } else NA_real_
            if (is.finite(pages) && pages > 0) return(pages * page_bytes / 1024^2)
            fallback_mb
        } else {
            fallback_mb
        }
    }, error = function(e) fallback_mb, warning = function(w) fallback_mb)
}

#' Draw `size` points from valid cells, density-weighted proportional to
#' `weight`'s own value (uniform over valid cells falls out automatically
#' when `weight` happens to be spatially constant - e.g. a plain 0/1
#' binary mask - via the `"degenerate"` branch below, not a separate calling
#' convention). Shared by `gm_convexity_index()`, `gm_span_index()`,
#' `gm_radial_concentration_index()` - the raster analogue of
#' `shapeindices::.sample_weighted_points()`. Every one of these indices is
#' a Monte Carlo estimate of an expectation over i.i.d. draws (mean
#' distance, mean fraction-outside, geometric median of an i.i.d. sample),
#' so this draw genuinely needs to BE i.i.d., not just "a
#' plausible-looking sample" - see each branch's own comments below for
#' how that's confirmed, not assumed, in each case.
#'
#' TWO DIFFERENT `terra::spatSample()` BACKENDS, chosen for a real memory
#' reason, confirmed by reading `terra`'s own source, not guessed:
#'  - `method = "random"` (used here whenever weight is spatially
#'    constant, and for each bin's own draw in the weighted/binned path
#'    below) bottoms out
#'    in `terra:::.sampleCellsRandom()`, which draws random CELL INDICES
#'    directly (`sample.int(ncell(x), ...)` - cheap integers, not
#'    materialised cell data) and reads only those specific cells' values
#'    via `terra:::add_cxyp()`'s `x[cnrs]` (targeted cell-number indexing,
#'    not a full-raster pull). Cost scales with `size` (times a rejection
#'    factor for `na.rm`), never with `n_valid_cells`.
#'  - `method = "weights"` bottoms out in `terra:::sampleWeights()`, which
#'    does `as.data.frame(x)` (materialising EVERY valid cell into an R
#'    data.frame - this is also what excludes invalid cells, since
#'    `as.data.frame()` drops `NA` rows) THEN
#'    `sample.int(nrow(res), size, prob = res[, ncol(res)], replace = replace)`.
#'    The `sample.int(..., replace = TRUE)` step is a perfectly good i.i.d.
#'    weighted-categorical sampler - the problem is the `as.data.frame()`
#'    step before it, whose cost scales with `n_valid_cells`, not `size`.
#'    Binning the WEIGHT VALUES first does NOT fix this on its own: the
#'    bottleneck is per-CELL materialisation, not the number of distinct
#'    values among those cells - a raster with a billion valid cells still
#'    needs a billion-row data.frame even if every one of those cells'
#'    weights has been rounded into one of 1000 bins.
#'
#' `replace = TRUE` throughout is load-bearing, not incidental, for BOTH
#' backends: `replace = FALSE` would make each successive draw depend on
#' every earlier one (sampling without replacement from a finite
#' population is NOT i.i.d.), silently breaking every index's own Monte
#' Carlo justification. Verified empirically (not just from source):
#' 50000 weighted draws from a small test raster reproduced the target
#' proportions to within 0.001, and the drawn sequence's own lag-1
#' autocorrelation was ~0 - consistent with genuine independence, not
#' merely "close to it in aggregate" (see `test-utils.R`'s own tests for
#' this, which run on every `R CMD check`, not just this one investigation).
#'
#' THE FIX for `method = "weights"`'s own memory scaling, in two parts:
#'  - a spatially-CONSTANT weight (the `"degenerate"` `.bin_index_raster()`
#'    mode - includes the common case of `weight` being a plain 0/1/NA
#'    binary mask) means every valid cell is equally likely regardless of
#'    `n_valid_cells`, so go straight to `method = "random"` - exact AND
#'    memory-safe at any raster size. This needs its own explicit branch,
#'    separate from `"exact"`: routing `"degenerate"` through
#'    `method = "weights"` instead would be fine when `n_valid_cells` is
#'    also small, but would silently reintroduce the FULL `as.data.frame()`
#'    cost for a huge raster that merely happens to have uniform weight -
#'    e.g. an ordinary large binary mask, a common case given there's no
#'    separate unweighted calling convention to reach for instead.
#'  - below `n_bins` valid cells (`"exact"` mode), `as.data.frame()`'s own
#'    cost is already bounded by `n_bins` regardless of total raster size,
#'    so `method = "weights"` there is cheap regardless.
#'  - otherwise (`"binned"` mode, a genuinely large raster with genuinely
#'    varying weight), `.bin_index_raster()` (same construction
#'    `.adaptive_density_bins()` uses for the weighted REFERENCE
#'    calculations, chunk-safe via `terra::app()` + `terra::freq()`)
#'    reduces the weight distribution to at most `n_bins` (density, count)
#'    strata, each treated as having constant density - genuinely the same
#'    piecewise-constant approximation the reference computations already
#'    accept, so no NEW kind of error is introduced. Each of the `size`
#'    draws is assigned a BIN i.i.d. from that small (`<= n_bins`-long)
#'    categorical distribution - trivial memory, `n_valid_cells` never
#'    enters into it - and then an ACTUAL cell is drawn from within that
#'    bin via `method = "random"` (real coordinates, not an approximated
#'    "bin-average" position), with a per-bin `exp` derived EXACTLY from
#'    that bin's own already-known cell count (a common bin needs almost
#'    no oversampling; a rare one gets proportionally more). See
#'    `.sample_from_bins()`'s own header for the residual edge case (a bin
#'    too rare for even a generous `exp` ceiling) and why its fallback is
#'    cheap precisely when it's needed.
#' @param valid logical terra SpatRaster (TRUE = inside the shape)
#' @param weight terra SpatRaster, the SAME raster the shape itself was
#'   derived from (`.valid_cells(weight)`) - need not already be masked to
#'   `valid`, masked internally
#' @param size integer, number of points to draw
#' @param n_bins integer, the exact/binned threshold used ONLY to bound
#'   `method = "weights"`'s own memory cost - an internal implementation
#'   detail, not exposed by any of the three index functions that call
#'   this (unlike the `n_bins` those functions expose for their own
#'   weighted REFERENCE computations, which trade off a different cost -
#'   see `gm_span_index()`'s own header for why that one is much smaller)
#' @return `size` x 2 numeric coordinate matrix, columns `x`, `y`
#' @noRd
.sample_valid_points <- function(valid, weight, size, n_bins = 1000) {
    b <- .bin_index_raster(valid, weight, n_bins)

    if (b$mode == "degenerate") {
        v <- terra::ifel(valid, 1, NA)
        s <- terra::spatSample(v, size = size, method = "random", xy = TRUE,
                                na.rm = TRUE, replace = TRUE, warn = FALSE)
        return(as.matrix(s[, c("x", "y")]))
    }
    if (b$mode == "exact") {
        w <- terra::ifel(valid, weight, NA)
        s <- terra::spatSample(w, size = size, method = "weights", xy = TRUE,
                                na.rm = TRUE, replace = TRUE, warn = FALSE)
        return(as.matrix(s[, c("x", "y")]))
    }

    .sample_from_bins(b, size)
}

#' Draw `size` points from a `"binned"` `.bin_index_raster()` result - one
#' BIN per point first (a cheap, `<= n_bins`-long categorical draw over
#' each bin's own total mass), then that point's actual coordinates via
#' `method = "random"` restricted to its own bin (see `.sample_valid_points()`'s
#' own header for why that backend, specifically, avoids the
#' `n_valid_cells`-scaling memory cost `method = "weights"` has). `exp`
#' (the rejection-sampling headroom `method = "random"` needs to reliably
#' find enough non-`NA` cells) is derived PER BIN from that bin's own
#' already-known cell count, not a fixed default - a common bin (most of
#' the raster) needs almost no oversampling; a rare one gets
#' proportionally more, capped at `exp_ceiling`. If a bin is rarer than
#' even that generous cap can reliably find, falls back to
#' `method = "weights"` restricted to JUST that bin - cheap precisely
#' because "too rare for rejection sampling to find easily" implies "few
#' actual member cells" for any realistic raster size, so materialising
#' only that bin's own (small) population is safe even though
#' materialising the WHOLE valid population would not be.
#'
#' ROW ORDER: results are assembled one bin at a time (all of one bin's
#' points, then the next bin's), NOT in the original per-point draw
#' order - and then explicitly shuffled before returning. This isn't
#' cosmetic: `gm_convexity_index()`/`gm_span_index()` pair up CONSECUTIVE rows
#' into (x1, x2), and real-world weight rasters routinely have spatial
#' autocorrelation (nearby cells tend to have similar density, hence land
#' in the same bin) - pairing same-bin rows without shuffling would
#' systematically pair spatially-nearby points together far more often
#' than genuine i.i.d. pairing would, biasing gm_span_index()'s own mean-
#' pairwise-distance estimate downward. Shuffling restores the
#' "any consecutive pair is a valid i.i.d. pair" property an
#' un-grouped i.i.d. sequence has for free - verified directly, not just
#' argued: `test-utils.R`'s own "consecutive pairs... show no spurious
#' within-pair correlation" test checks `corr(x1_x, x2_x)` stays ~0 even
#' under a strongly spatially-correlated weight (a left-to-right
#' gradient) and aggressively coarse binning - the scenario where a
#' broken shuffle would show up as a spurious positive correlation, not
#' as "close to independent".
#' @param b a `"binned"`-mode `.bin_index_raster()` result
#' @param size number of points to draw
#' @param safety_margin how many EXPECTED hits per bin to aim for (not
#'   just 1) before considering that bin's own rejection sampling reliable
#' @param exp_ceiling upper bound on the per-bin oversampling factor,
#'   above which the `method = "weights"` fallback kicks in instead
#' @return `size` x 2 numeric coordinate matrix, columns `x`, `y`
#' @noRd
.sample_from_bins <- function(b, size, safety_margin = 3, exp_ceiling = 2000) {
    f <- terra::freq(b$r_idx)
    bin_value <- f$value
    bin_count <- f$count
    bin_mass <- (b$rng[1] + (bin_value + 0.5) * b$bin_width) * bin_count
    n_bin_actual <- length(bin_value)
    n_grid_cells <- terra::ncell(b$r_idx)

    bin_draw <- sample.int(n_bin_actual, size, replace = TRUE, prob = bin_mass)
    needed <- tabulate(bin_draw, nbins = n_bin_actual)

    # per-bin rejection sampling, done DIRECTLY here (sample.int() over
    # raw cell numbers + a targeted b$r_idx[cnrs] read) rather than via
    # terra::spatSample(method = "random") - that path was tried first and
    # hits a real terra bug (terra:::.sampleCellsMemory()'s `v[i, 1]`
    # silently drops to a bare vector when exactly one cell matches,
    # crashing sample.int() with "length(n) == 1L is not TRUE") on small-
    # to-moderate rasters with a sparsely-populated bin - reproduced and
    # confirmed via a direct call before switching to this lower-level
    # equivalent, which sidesteps that code path entirely
    collected <- vector("list", n_bin_actual)
    for (k in which(needed > 0)) {
        m_k <- needed[k]
        frac_k <- bin_count[k] / n_grid_cells
        exp_k <- min(exp_ceiling, max(5, ceiling(safety_margin / frac_k)))

        hits <- integer(0)
        exp_try <- exp_k
        for (round in 1:5) {
            cnrs <- sample.int(n_grid_cells, m_k * exp_try, replace = TRUE)
            vals <- b$r_idx[cnrs][[1]]
            hits <- c(hits, cnrs[!is.na(vals) & vals == bin_value[k]])
            if (length(hits) >= m_k) break
            exp_try <- exp_try * 4
        }
        if (length(hits) < m_k) {
            # still short after widening exp_try 4x for 5 rounds - this
            # bin is rarer than even a generous budget can reliably find;
            # materialise JUST this bin directly instead (cheap precisely
            # because "hard to find at random" implies "few actual member
            # cells", for any raster size)
            bin_mask <- terra::ifel(b$r_idx == bin_value[k], 1, NA)
            top_up <- terra::spatSample(bin_mask, size = m_k - length(hits), method = "weights",
                                         xy = TRUE, na.rm = TRUE, replace = TRUE, warn = FALSE)
            hits <- c(hits, terra::cellFromXY(b$r_idx, as.matrix(top_up[, c("x", "y")])))
        }
        collected[[k]] <- hits[seq_len(m_k)]
    }

    all_cells <- unlist(collected)
    xy <- terra::xyFromCell(b$r_idx, all_cells)
    xy[sample.int(nrow(xy)), , drop = FALSE]
}

#' The largest Monte Carlo `size` considered safe to sample, given a
#' memory budget - shared by `gm_convexity_index()`/`gm_span_index()`/
#' `gm_radial_concentration_index()`. Mirrors
#' `shapeindices::.safe_deterministic_tri_ceiling()`'s "invert the actual
#' cost formula against `.available_memory_mb()`" approach, adapted to
#' Monte Carlo sample size instead of triangle count.
#' @param valid logical terra SpatRaster - only its own dimensions matter
#'   here, not its actual TRUE/FALSE pattern
#' @param formula which index's own cost this ceiling protects: `"point"`
#'   (O(size) - `gm_span_index()`'s/`gm_radial_concentration_index()`'s own
#'   sampled point cloud), `"line"` (O(size * raster diagonal in cells) -
#'   `gm_convexity_index()`'s own per-line raster discretisation, the more
#'   expensive of the two, since each line needs its own raster-resolution
#'   sub-sampling, not just its two endpoints), or `"geodesic"` (see below)
#' @param bytes_per_unit conservative bytes needed per point (`formula =
#'   "point"`) or per discretised line sub-point (`formula = "line"`) -
#'   unused for `formula = "geodesic"`, which has its own, differently-
#'   scaled constant (see below)
#' @param mem_fraction fraction of estimated available memory to budget
#' @return integer >= 1, the largest `size` (or, for `"geodesic"`, `K`)
#'   considered safe
#' @noRd
.safe_mc_size_ceiling <- function(valid, formula = c("point", "line", "geodesic"), bytes_per_unit = 64, mem_fraction = 0.2) {
    formula <- match.arg(formula)
    budget_bytes <- .available_memory_mb() * 1024^2 * mem_fraction
    max_units <- budget_bytes / bytes_per_unit
    if (formula == "point") return(max(1L, as.integer(floor(max_units))))

    if (formula == "geodesic") {
        # gm_geodesic_span_index()/gm_geodesic_chord_index()'s own K
        # sampled points each cost one WHOLE-RASTER terra::gridDist()
        # call (K calls total, run sequentially, each result discarded
        # once its own K-1 pairwise distances are read off) - genuinely a
        # different cost shape from "point"/"line" above: each individual
        # call is itself already a chunk-safe, terra-native operation (the
        # same way gm_depth_index()'s own single terra::distance() call
        # needs no ceiling at all), so the real resource being protected
        # here is WALL-CLOCK TIME (K sequential whole-raster sweeps), not
        # peak memory. Still expressed as a memory-budget-shaped
        # ceiling, reusing `.available_memory_mb()` only as this
        # package's own already-established proxy for "how much
        # computation this machine can reasonably absorb" - not a literal
        # per-cell memory accounting the way "point"/"line" are. The
        # constant below (not `bytes_per_unit`, which is scaled for a
        # genuinely different unit) is calibrated directly against a
        # measured benchmark (~5e-8 seconds/cell/call on ordinary
        # hardware - K=100 on a 1.44M-cell raster took ~7s), targeting a
        # ceiling in the single-digit-to-low-tens-of-seconds range at the
        # `mem_fraction`-scaled budget rather than a true memory bound.
        n_cells <- terra::ncell(valid)
        geodesic_unit_bytes <- 2
        max_geo_units <- budget_bytes / geodesic_unit_bytes
        return(max(1L, as.integer(floor(max_geo_units / n_cells))))
    }

    # "line": each of size/2 lines needs ~diag_cells sub-points at raster
    # resolution (half a cell diagonal step) - diag_cells (the raster's
    # own extent) is a conservative upper bound on any sampled line's
    # length, since every sampled point lies within that extent
    diag_cells <- sqrt(sum(dim(valid)[1:2]^2))
    max(1L, as.integer(floor(2 * max_units / diag_cells)))
}

#' Hard-stop (not a silent clamp) if `size` exceeds the memory-derived
#' ceiling - consistent with every other ceiling in this package.
#' @param arg_name the caller's own argument name for `size` in its error
#'   message - `"size"` for every caller except `gm_geodesic_span_index()`/
#'   `gm_geodesic_chord_index()`, which call this with their own
#'   `n_points` value under a deliberately different argument name (see
#'   R/geodesic-index.R's own file header for why) and need the error
#'   text to say so, not "size"
#' @noRd
.check_mc_size <- function(size, valid, formula, fn_name, arg_name = "size") {
    ceiling_size <- .safe_mc_size_ceiling(valid, formula = formula)
    if (size > ceiling_size) {
        stop(fn_name, "(): `", arg_name, "` (", size, ") exceeds the estimated safe ceiling (",
             ceiling_size, ") given available memory and this raster's extent. ",
             "Lower `", arg_name, "`, or free up memory before calling again - not silently clamped.")
    }
}

#' Fraction of each line's length lying outside `valid`, via raster
#' discretisation - the raster analogue of `shapeindices`' vector
#' segment-vs-boundary-edge line clipping, but simpler: no vector geometry
#' at all, just sample each line at raster resolution and look up whether
#' each sub-point falls on a valid cell. Every line gets its OWN step
#' count derived from its own physical length and the raster's own cell
#' size (finer sampling buys nothing past the mask's own resolution),
#' then all lines' sub-points are flattened into one big coordinate
#' matrix so `terra::extract()` runs once, not once per line.
#' @param valid logical terra SpatRaster (TRUE = inside the shape)
#' @param x1,x2 Nx2 coordinate matrices, line endpoints (line `i` runs
#'   from `x1[i, ]` to `x2[i, ]`)
#' @return numeric vector, length N, fraction of each line's length where
#'   the raster is NOT valid (`NA` counts as outside, same as a hole)
#' @noRd
.frac_outside_line <- function(valid, x1, x2) {
    n <- nrow(x1)
    len <- sqrt((x2[, 1] - x1[, 1])^2 + (x2[, 2] - x1[, 2])^2)
    step <- min(terra::res(valid)) / 2
    n_sub <- pmax(2L, ceiling(len / step) + 1L)

    line_id <- rep(seq_len(n), n_sub)
    tt <- unlist(lapply(n_sub, function(k) seq(0, 1, length.out = k)))
    px <- x1[line_id, 1] + tt * (x2[line_id, 1] - x1[line_id, 1])
    py <- x1[line_id, 2] + tt * (x2[line_id, 2] - x1[line_id, 2])

    inside <- terra::extract(valid, cbind(px, py))[[1]]
    inside[is.na(inside)] <- FALSE

    # uniform t-spacing along each line => the plain mean over that line's
    # own sub-points already IS its length-weighted mean; no extra
    # per-point weighting needed
    as.numeric(tapply(1 - as.numeric(inside), line_id, mean))
}

#' Concentric-rings D1 reference (distance-to-median / distance-to-
#' geometric-median-equivalent), raster version - each ring is either one
#' valid cell (exact path) or one density bin spanning several cells
#' (binned path, see `.adaptive_density_bins()` above), so ring k's outer
#' boundary comes from *cumulative cell count* rather than a fixed
#' `k * cell_area`. Shared by `gm_depth_index()` (as `R - D1_ref`) and
#' `gm_radial_concentration_index()` (directly, as `D1_ref` itself) - the
#' same closed-form mean-radius-to-centre integral either way, exact
#' (polynomial, not quadrature) regardless of ring count, unlike
#' `gm_span_index()`'s own `|x-y|`-kernel reference.
#' @param cell_area scalar, the (identical) area of one raster cell
#' @param density_sorted numeric vector, density per ring, sorted
#'   descending
#' @param count numeric vector, number of cells in each ring (default: one
#'   cell per ring, i.e. `density_sorted` itself is the per-cell vector)
#' @return D1_ref
#' @noRd
.annulus_reference_D1_raster <- function(cell_area, density_sorted, count = rep(1, length(density_sorted))) {
    S <- cell_area * cumsum(count)
    S_prev <- c(0, S[-length(S)])
    r_hi <- sqrt(S / pi)
    r_lo <- sqrt(S_prev / pi)
    mass <- density_sorted * count
    W <- mass / sum(mass)
    mean_r <- (2 / 3) * (r_lo^2 + r_lo * r_hi + r_hi^2) / (r_lo + r_hi)
    sum(W * mean_r)
}

#' Cleared on exit of every exported index function via
#' `on.exit(.cleanup_tmpfiles(), add = TRUE)` - terra spills intermediate
#' rasters to disk under memory pressure and does not always clean up
#' eagerly mid-session (only guaranteed at normal R session end), which can
#' genuinely fill a disk for a package meant to run in loops over many
#' rasters. `.onLoad()` points `terraOptions(tempdir=)` at a dedicated
#' gridmorph subdirectory so this cleanup never sweeps up unrelated
#' concurrent terra usage elsewhere in the same R session.
#' @noRd
.cleanup_tmpfiles <- function() {
    terra::tmpFiles(remove = TRUE)
}
