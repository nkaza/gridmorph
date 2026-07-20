## =========================================================================
## Running all indices (or a chosen subset) over one raster in one call -
## the raster analogue of shapeindices::shape_indices(). Same `which`-
## driven design, same short index names (so results from the two
## packages line up directly for comparison), same named-numeric-vector
## return shape.
##
## NO SHARED-MESH OPTIMISATION - unlike shapeindices' own wrapper (which
## triangulates once and reuses that CDT mesh across every requested
## index, a genuinely expensive shared cost), gridmorph has no analogous
## expensive shared setup: `.valid_cells()`/`.mass_raster()` are cheap
## per-call raster derivations, not worth hoisting out. What IS worth
## doing once, not once per requested index: the planar-CRS check -
## calling all 13 index functions independently on a missing-CRS raster
## would otherwise print the same "no CRS" warning 13 times.

## canonical order - also what "all" expands to, and what the returned
## vector's own names appear in regardless of the order `which` was given
.ALL_GM_INDICES <- c("depth", "moment_of_inertia", "moment_isotropy", "directional_balance",
                     "convexity", "span", "radial_concentration", "hull_ratio", "polsby_popper",
                     "width_length_ratio", "reock", "detour", "exchange")

## the six classic metrics take no weighted/size/seed/n_bins arguments at
## all - see classical-metrics.R's own file header for why
.GM_CLASSIC_INDICES <- c("hull_ratio", "polsby_popper", "width_length_ratio", "reock", "detour", "exchange")

#' Filters `dots` down to the names a function actually accepts, so one
#' shared `...` can feed all thirteen index functions with only partly-
#' overlapping arguments (e.g. `weighted`/`n_bins`, not accepted by the
#' six classic metrics at all) without `do.call()` erroring on an unused
#' argument.
#' @param fn the target function
#' @param dots a named list
#' @return the subset of `dots` whose names are in `names(formals(fn))`
#' @noRd
.dots_for <- function(fn, dots) dots[names(dots) %in% names(formals(fn))]

#' Resolves `gm_shape_indices()`'s own `which` argument.
#' @param which `"all"`, or a character vector naming a subset of
#'   `.ALL_GM_INDICES`
#' @return character vector, the requested subset in canonical order
#' @noRd
.resolve_which <- function(which) {
    if (identical(which, "all")) return(.ALL_GM_INDICES)
    if (!is.character(which)) {
        stop("`which` must be \"all\" or a character vector of index names; got ",
             class(which)[1], ".")
    }
    unknown <- setdiff(which, .ALL_GM_INDICES)
    if (length(unknown) > 0) {
        stop("Unknown index name(s) in `which`: ", paste(unknown, collapse = ", "),
             ". Valid choices: ", paste(.ALL_GM_INDICES, collapse = ", "), ", or \"all\".")
    }
    .ALL_GM_INDICES[.ALL_GM_INDICES %in% which]
}

#' All indices (or a chosen subset) for a single terra raster
#'
#' @param rast a terra SpatRaster - see [gm_depth_index()] for the shared
#'   shape/hole conventions every index in this package uses.
#' @param which `"all"` (default), or a character vector naming a subset
#'   of these thirteen values - each listed here with the function it
#'   actually calls, since the `which` string and the function name aren't
#'   identical:
#'
#'   * `"depth"` - [gm_depth_index()]
#'   * `"moment_of_inertia"` - [gm_moment_of_inertia_index()]
#'   * `"moment_isotropy"` - [gm_moment_isotropy_index()]
#'   * `"directional_balance"` - [gm_directional_balance_index()]
#'   * `"convexity"` - [gm_convexity_index()]
#'   * `"span"` - [gm_span_index()]
#'   * `"radial_concentration"` - [gm_radial_concentration_index()]
#'   * `"hull_ratio"` - [gm_hull_ratio_index()]
#'   * `"polsby_popper"` - [gm_polsby_popper_index()]
#'   * `"width_length_ratio"` - [gm_width_length_ratio_index()]
#'   * `"reock"` - [gm_reock_index()]
#'   * `"detour"` - [gm_detour_index()]
#'   * `"exchange"` - [gm_exchange_index()]
#'
#'   The same short names `shapeindices::shape_indices()` uses, so results
#'   from the two packages line up directly. An unrecognised name in
#'   `which` errors immediately, listing all thirteen valid values.
#' @param ... passed to whichever of the thirteen index functions actually
#'   accept each named argument (e.g. `weighted` is accepted by the first
#'   seven, silently ignored for the six classic metrics, which have no
#'   weighted form at all - see their own file header for why; `size`/
#'   `seed` are accepted only by the three Monte Carlo indices
#'   (`gm_convexity_index()`/`gm_span_index()`/`gm_radial_concentration_index()`);
#'   `n_bins` is accepted by `gm_depth_index()`/
#'   `gm_moment_of_inertia_index()`/`gm_span_index()`/
#'   `gm_radial_concentration_index()` - passing it explicitly overrides
#'   each of their own individual defaults, including `gm_span_index()`'s
#'   deliberately smaller one).
#' @return named numeric vector, one entry per requested index, in
#'   canonical order
#' @examples
#' r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
#' terra::values(r) <- 0
#' r[10:30, 10:30] <- 1
#' gm_shape_indices(r, size = 1000, seed = 1)
#'
#' # a subset
#' gm_shape_indices(r, which = c("hull_ratio", "polsby_popper", "reock"))
#' @export
gm_shape_indices <- function(rast, which = "all", ...) {
    which <- .resolve_which(which)
    .check_planar_crs(rast, "gm_shape_indices")
    dots <- list(...)
    # weighted defaults to TRUE, matching every weighted-capable index's
    # own default - warn ONCE here, up front, same treatment as the
    # planar-CRS check just above, rather than once per weighted-capable
    # index called below
    weighted <- if (is.null(dots$weighted)) TRUE else dots$weighted
    if (weighted && terra::is.factor(rast)) .warn_categorical_weighted(rast)

    out <- c()
    # every individual gm_*_index() call re-checks the planar CRS (and,
    # for the seven weighted-capable ones, the categorical-raster
    # condition) on its own - already told the caller once, above; muffle
    # the redundant repeats of both (missing-CRS warnings and the
    # categorical-weighted warning only - a geographic CRS already
    # stopped everything at the .check_planar_crs() call above)
    .muffle_crs_warnings(.muffle_categorical_warnings({
        if ("depth" %in% which) {
            out["depth"] <- do.call(gm_depth_index, c(list(rast = rast), .dots_for(gm_depth_index, dots)))$index
        }
        if ("moment_of_inertia" %in% which) {
            out["moment_of_inertia"] <- do.call(gm_moment_of_inertia_index, c(list(rast = rast), .dots_for(gm_moment_of_inertia_index, dots)))$index
        }
        if ("moment_isotropy" %in% which) {
            out["moment_isotropy"] <- do.call(gm_moment_isotropy_index, c(list(rast = rast), .dots_for(gm_moment_isotropy_index, dots)))$index
        }
        if ("directional_balance" %in% which) {
            out["directional_balance"] <- do.call(gm_directional_balance_index, c(list(rast = rast), .dots_for(gm_directional_balance_index, dots)))$index
        }
        if ("convexity" %in% which) {
            out["convexity"] <- do.call(gm_convexity_index, c(list(rast = rast), .dots_for(gm_convexity_index, dots)))$index
        }
        if ("span" %in% which) {
            out["span"] <- do.call(gm_span_index, c(list(rast = rast), .dots_for(gm_span_index, dots)))$index
        }
        if ("radial_concentration" %in% which) {
            out["radial_concentration"] <- do.call(gm_radial_concentration_index, c(list(rast = rast), .dots_for(gm_radial_concentration_index, dots)))$index
        }
        if ("hull_ratio" %in% which) out["hull_ratio"] <- gm_hull_ratio_index(rast)$index
        if ("polsby_popper" %in% which) out["polsby_popper"] <- gm_polsby_popper_index(rast)$index
        if ("width_length_ratio" %in% which) out["width_length_ratio"] <- gm_width_length_ratio_index(rast)$index
        if ("reock" %in% which) out["reock"] <- gm_reock_index(rast)$index
        if ("detour" %in% which) out["detour"] <- gm_detour_index(rast)$index
        if ("exchange" %in% which) out["exchange"] <- gm_exchange_index(rast)$index
    }))

    out[which]
}
