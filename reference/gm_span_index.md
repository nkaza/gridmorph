# Span index for a terra raster

The raster analogue of
[`shapeindices::span_index()`](https://nkaza.github.io/shapeindices/reference/span_index.html).
`D_ref/D`, in `(0, 1]`, where `D` is the mean distance between two
random interior points and `D_ref` is the same quantity for the
reference shape (a circle, unweighted; a concentric annulus, weighted) -
both provable minimisers of `D`. Distinct from
[`gm_moment_of_inertia_index()`](https://nkaza.github.io/gridmorph/reference/gm_moment_of_inertia_index.md),
which a squared-distance version of this index would just collapse to.

## Usage

``` r
gm_span_index(rast, weighted = TRUE, size = 3000, seed = NULL, n_bins = 100)
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

  number of points to sample (paired up consecutively into `size %/% 2`
  point pairs) - matches
  [`terra::spatSample()`](https://rspatial.github.io/terra/reference/sample.html)'s
  own argument name and meaning directly. Checked against a
  memory-derived ceiling before running; hard-stops, not a silent clamp,
  if exceeded.

- seed:

  optional RNG seed

- n_bins:

  integer, the exact/binned threshold and bin count for the weighted
  reference's concentric-rings construction. Default `100`, much lower
  than
  [`gm_depth_index()`](https://nkaza.github.io/gridmorph/reference/gm_depth_index.md)'s/[`gm_moment_of_inertia_index()`](https://nkaza.github.io/gridmorph/reference/gm_moment_of_inertia_index.md)'s
  `1000` - see file header for why (this reference's own cost is
  quadratic in ring count, not linear).

## Value

`list(index, D, D_ref, area, n_valid_cells)`

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_span_index(r, size = 2000, seed = 1)$index
#> [1] 0.9833546
```
