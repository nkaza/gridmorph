# gridmorph: Raster-Native Shape Indices and Morphological Operators

Shape/compactness/dispersal indices for a single terra SpatRaster,
treating a whole shape - however many disjoint patches or holes it has -
as one combined shape, the same way `shapeindices` treats a
MULTIPOLYGON. The shape itself is derived directly from the raster's own
values (non-`NA`, non-zero cells); a `weighted` argument on every index
controls whether those values are also used as a density/mass field
(`TRUE`, default) or ignored beyond defining the shape (`FALSE`). Also
exposes morphological operators (erosion, dilation, opening, closing,
top-hat, bottom-hat) directly on `terra` SpatRaster objects.

## Planar rasters only

Every index in this package assumes `rast` uses PLANAR (projected)
coordinates. This is *not* because area, distance, or perimeter can't be
computed correctly on geographic (longitude/latitude) data in general -
`terra`'s own
[`terra::expanse()`](https://rspatial.github.io/terra/reference/expanse.html),
[`terra::perim()`](https://rspatial.github.io/terra/reference/perim.html),
and
[`terra::distance()`](https://rspatial.github.io/terra/reference/distance.html)
all do that correctly, via geodesic math, even in a geographic CRS (see
`terra`'s own documentation, which explicitly notes that area is often
*better* computed directly in longitude/latitude than by reprojecting
first). The real problem is narrower: this package's own point-based
geometry - line-crossing tests for
[`gm_convexity_index()`](https://nkaza.github.io/gridmorph/reference/gm_convexity_index.md),
pairwise-distance accumulation for
[`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)/
[`gm_radial_concentration_index()`](https://nkaza.github.io/gridmorph/reference/gm_radial_concentration_index.md),
moment-tensor calculations for
[`gm_moment_of_inertia_index()`](https://nkaza.github.io/gridmorph/reference/gm_moment_of_inertia_index.md)
and its relatives, and minimum-enclosing- circle fitting for
[`gm_reock_index()`](https://nkaza.github.io/gridmorph/reference/gm_reock_index.md)/[`gm_detour_index()`](https://nkaza.github.io/gridmorph/reference/gm_detour_index.md) -
works directly in raw x/y coordinate space throughout, treating it as
Cartesian, and is not yet geodesic-aware. Every index function checks
this and **errors** if `rast` has a real geographic CRS - reproject to a
suitable projected CRS first (e.g. via
[`terra::project()`](https://rspatial.github.io/terra/reference/project.html)).
A raster with no CRS at all is allowed (warned about, not blocked) and
treated as already planar, matching the common case of a purely
synthetic or abstract grid with no real-world geographic meaning. The
morphological operators
([`erode()`](https://nkaza.github.io/gridmorph/reference/erode.md),
[`dilate()`](https://nkaza.github.io/gridmorph/reference/dilate.md),
[`opening()`](https://nkaza.github.io/gridmorph/reference/opening.md),
[`closing()`](https://nkaza.github.io/gridmorph/reference/closing.md),
[`tophat()`](https://nkaza.github.io/gridmorph/reference/tophat.md),
[`bottomhat()`](https://nkaza.github.io/gridmorph/reference/bottomhat.md))
are the one exception - they only reason about cell adjacency via a
structuring-element kernel, never physical distance, so they work
identically regardless of CRS.

## See also

Useful links:

- <https://github.com/nkaza/gridmorph>

- Report bugs at <https://github.com/nkaza/gridmorph/issues>

## Author

**Maintainer**: Nikhil Kaza <kaza@cs.unc.edu>

Authors:

- Nikhil Kaza <kaza@cs.unc.edu>
