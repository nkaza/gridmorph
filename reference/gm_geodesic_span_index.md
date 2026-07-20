# Geodesic span index for a terra raster

Mean GEODESIC distance between two random interior points, relative to
an equal-area circle - the raster analogue of a shape-confined-path
variant of
[`shapeindices::span_index()`](https://nkaza.github.io/shapeindices/reference/span_index.html)
(see file header; deferred in the vector package for lacking cheap
shortest-path machinery, not needed here). `index = D_ref/D`, in
`(0, 1]`, `= 1` iff the shape is (almost everywhere) a disk. Distinct
from
[`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)
(same construction, Euclidean straight-line distance instead) - the two
agree exactly on any CONVEX shape (where the straight line between two
interior points never leaves it) and diverge only once concavity forces
a detour.

## Usage

``` r
gm_geodesic_span_index(
  rast,
  weighted = TRUE,
  n_points = 40,
  seed = NULL,
  n_bins = 100
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

- n_points:

  number of points to sample (ALL pairs among them are used - see above,
  this is NOT
  [`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)'s
  own `size` convention, deliberately a different argument name so the
  two never collide through
  [`gm_shape_indices()`](https://nkaza.github.io/gridmorph/reference/gm_shape_indices.md)'s
  own shared `...`). Checked against a memory/time-derived ceiling
  before running.

- seed:

  optional RNG seed

- n_bins:

  integer, the exact/binned threshold and bin count for the weighted
  reference's concentric-rings construction - see
  [`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)'s
  own doc for the same parameter's meaning

## Value

`list(index, D, D_ref, area, n_valid_cells)`

## Details

`n_points` means something DIFFERENT here than `size` does in
[`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md) -
a DELIBERATELY DIFFERENT ARGUMENT NAME, not an oversight, precisely to
avoid this - see file header. `n_points = K` draws K points and averages
ALL K(K-1)/2 pairwise distances among them (not K/2 independent pairs),
so a much smaller `n_points` already gives a comparable number of
pairs - but each point costs one whole-raster
[`terra::gridDist()`](https://rspatial.github.io/terra/reference/gridDist.html)
call, a substantially higher per-point cost than
[`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)'s
own closed-form Euclidean distance. Checked against a
memory/time-derived ceiling before running (`formula = "geodesic"`,
R/utils.R); hard-stops, not a silent clamp, if exceeded.

BOTH `D` and `D_ref` carry real Monte Carlo noise - more than
[`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)'s
own estimate at a comparable sample size, since
[`terra::gridDist()`](https://rspatial.github.io/terra/reference/gridDist.html)'s
own angular quantization adds variability on top of ordinary
point-sampling noise (verified: several percent spread across seeds even
at `n_points = 200` on a test disk). Increase `n_points` for a more
precise answer; there is no way to eliminate this noise entirely the way
the closed-form Euclidean reference has none.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_geodesic_span_index(r, n_points = 25, seed = 1)$index
#> [1] 0.9199875
```
