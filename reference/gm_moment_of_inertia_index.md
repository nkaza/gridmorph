# Polar moment-of-inertia compactness/dispersal index for a terra raster

The raster analogue of
[`shapeindices::moment_of_inertia_index()`](https://nkaza.github.io/shapeindices/reference/moment_of_inertia_index.html).
index in `(0, 1]`, `= 1` iff the shape is (almost everywhere) a disk.

## Usage

``` r
gm_moment_of_inertia_index(rast, weighted = TRUE, n_bins = 1000)
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

- n_bins:

  integer, the exact/binned threshold and bin count for the weighted
  reference's concentric-rings construction (only matters when
  `weighted = TRUE` and the density genuinely varies -
  `weighted = FALSE` always collapses to a single ring regardless of
  `n_bins`) - exact (no accuracy cost) when the number of valid cells is
  `<= n_bins`, a K-bin histogram approximation above that. Default
  `1000` measured at ~0.1% relative error on a genuinely non-uniform
  weight pattern; see `.adaptive_density_bins()` in `R/utils.R`.

## Value

list(index, J, Ixx, Iyy, Ixy, J_ref, area, centroid)
