#' gridmorph: Raster-Native Shape Indices and Morphological Operators
#'
#' Shape/compactness/dispersal indices for a single terra SpatRaster,
#' treating a whole shape - however many disjoint patches or holes it has -
#' as one combined shape, the same way `shapeindices` treats a
#' MULTIPOLYGON. The shape itself is derived directly from the raster's own
#' values (non-`NA`, non-zero cells); a `weighted` argument on every index
#' controls whether those values are also used as a density/mass field
#' (`TRUE`, default) or ignored beyond defining the shape (`FALSE`). Also
#' exposes morphological operators (erosion, dilation, opening, closing,
#' top-hat, bottom-hat) directly on `terra` SpatRaster objects.
#'
#' # Planar rasters only
#'
#' Every index in this package assumes `rast` uses PLANAR (projected)
#' coordinates. This is *not* because area, distance, or perimeter can't be
#' computed correctly on geographic (longitude/latitude) data in
#' general - `terra`'s own [terra::expanse()], [terra::perim()], and
#' [terra::distance()] all do that correctly, via geodesic math, even in a
#' geographic CRS (see `terra`'s own documentation, which explicitly notes
#' that area is often *better* computed directly in longitude/latitude than
#' by reprojecting first). The real problem is narrower: this package's own
#' point-based geometry - line-crossing tests for `gm_convexity_index()`,
#' pairwise-distance accumulation for `gm_span_index()`/
#' `gm_radial_concentration_index()`, moment-tensor calculations for
#' `gm_moment_of_inertia_index()` and its relatives, and minimum-enclosing-
#' circle fitting for `gm_reock_index()`/`gm_detour_index()` - works
#' directly in raw x/y coordinate space throughout, treating it as
#' Cartesian, and is not yet geodesic-aware. Every index function checks
#' this and **errors** if `rast` has a real geographic CRS - reproject to a
#' suitable projected CRS first (e.g. via [terra::project()]). A raster
#' with no CRS at all is allowed (warned about, not blocked) and treated as
#' already planar, matching the common case of a purely synthetic or
#' abstract grid with no real-world geographic meaning. The morphological
#' operators (`erode()`, `dilate()`, `opening()`, `closing()`, `tophat()`,
#' `bottomhat()`) are the one exception - they only reason about cell
#' adjacency via a structuring-element kernel, never physical distance, so
#' they work identically regardless of CRS.
#'
#' @keywords internal
#' @import terra
"_PACKAGE"

#' @noRd
.onLoad <- function(libname, pkgname) {
    terra::terraOptions(parallel = TRUE)
    gm_tempdir <- file.path(tempdir(), "gridmorph")
    dir.create(gm_tempdir, showWarnings = FALSE, recursive = TRUE)
    terra::terraOptions(tempdir = gm_tempdir)
}
