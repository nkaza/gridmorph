# Convexity/dispersal index of a terra raster

The raster analogue of
[`shapeindices::convexity_index()`](https://nkaza.github.io/shapeindices/reference/convexity_index.html).
1 minus the expected fraction of a random interior line lying outside
the shape. index in `[0, 1]`, `1` = convex; lower means more concave
and/or more spatially dispersed. Handles holes and multi-part shapes
with no special-casing: a line between two disjoint parts is mostly
"outside", exactly like a vector MULTIPOLYGON's own line-clipping.

## Usage

``` r
gm_convexity_index(rast, weighted = TRUE, size = 3000, seed = NULL)
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
  lines) - matches
  [`terra::spatSample()`](https://rspatial.github.io/terra/reference/sample.html)'s
  own argument name and meaning directly (a raw point count, not a line
  count). Checked against a memory-derived ceiling before running (see
  `.safe_mc_size_ceiling()` in `R/utils.R`); hard-stops, not a silent
  clamp, if exceeded.

- seed:

  optional RNG seed

## Value

`list(index, n_lines, area, n_valid_cells)`

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_convexity_index(r, size = 500, seed = 1)$index
#> [1] 1
```
