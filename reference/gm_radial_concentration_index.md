# Radial concentration index for a terra raster

The raster analogue of
[`shapeindices::radial_concentration_index()`](https://nkaza.github.io/shapeindices/reference/radial_concentration_index.html).
`D1_ref/D1`, in `(0, 1]`, where `D1` is the mean distance from random
interior points to the shape's own geometric median (the point
minimising that mean distance - not the centroid) and `D1_ref` is the
same quantity for the reference shape (a circle, unweighted; a
concentric annulus, weighted) - both provable minimisers of `D1`.

## Usage

``` r
gm_radial_concentration_index(
  rast,
  weighted = TRUE,
  size = 3000,
  seed = NULL,
  n_bins = 1000
)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

- weighted:

  logical, default `TRUE`. When `TRUE`, `rast`'s own cell values are
  used as the density/mass throughout. When `FALSE`, every valid cell is
  treated as equally massed regardless of `rast`'s actual values -
  exactly reproducing the plain (unweighted) index even if `rast` is
  itself a continuous raster. If `rast` is a CATEGORICAL raster
  (`terra::is.factor(rast)` is `TRUE` - e.g. a land-use classification
  with an attached levels table) and `weighted = TRUE`, the category
  CODE NUMBERS themselves get used as mass, which is rarely the intended
  comparison (a class coded `3` outweighs one coded `1` for no
  meaningful reason) - a warning is issued but the computation still
  proceeds on the literal codes, since `weighted = TRUE` was explicitly
  requested; pass `weighted = FALSE` for categorical rasters unless the
  codes really are meant as magnitudes. This check never fires on
  non-categorical (including plain binary 0/1) rasters, so
  `weighted = TRUE` is always safe there.

- size:

  number of points to sample and run Weiszfeld's algorithm on directly
  (no pairing, unlike
  [`gm_convexity_index()`](https://nkaza.github.io/gridmorph/reference/gm_convexity_index.md)/[`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)) -
  matches
  [`terra::spatSample()`](https://rspatial.github.io/terra/reference/sample.html)'s
  own argument name and meaning directly. Checked against a
  memory-derived ceiling before running; hard-stops, not a silent clamp,
  if exceeded.

- seed:

  optional RNG seed

- n_bins:

  integer, the exact/binned threshold and bin count for the weighted
  reference's concentric-rings construction. Default `1000`, same as
  [`gm_depth_index()`](https://nkaza.github.io/gridmorph/reference/gm_depth_index.md)'s -
  this reference is an exact closed-form sum, not quadrature, so (unlike
  [`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md))
  there's no O(rings^2) cost to guard against.

## Value

`list(index, D1, D1_ref, area, center, n_valid_cells)`. `center` is the
geometric median found by Weiszfeld's algorithm - may be non-unique for
a symmetric multi-part shape (see `.geometric_median()`'s own comments),
in which case it can land anywhere along the minimising segment,
including inside a hole or the gap between multi-part pieces; the index
value itself is unaffected.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_radial_concentration_index(r, size = 2000, seed = 1)$index
#> [1] 0.9954244
```
