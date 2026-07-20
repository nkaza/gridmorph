# Geodesic chord index for a terra raster

Mean GEODESIC distance between two random points on the shape's own
BOUNDARY, relative to an equal-area circle - the raster analogue of
Angel, Parent & Civco (2010)'s "Traversal Index" (see file header for
the naming choice and why it's called "chord" here, not "traversal").
`index = D_ref/D`, in `(0, 1]`, `= 1` iff the shape is (almost
everywhere) a disk.

## Usage

``` r
gm_geodesic_chord_index(rast, n_points = 40, seed = NULL)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

- n_points:

  number of boundary points to sample (ALL pairs among them are used -
  see
  [`gm_geodesic_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_span_index.md)'s
  own doc for why this argument is named differently from
  [`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)'s
  `size`). Checked against a memory/time-derived ceiling before running.

- seed:

  optional RNG seed

## Value

`list(index, D, D_ref, area, n_valid_cells, n_boundary_cells)`

## Details

No `weighted` argument - boundary cells carry no interior mass to weight
by (see file header), matching
[`gm_hull_ratio_index()`](https://nkaza.github.io/gridmorph/reference/gm_hull_ratio_index.md)/
[`gm_polsby_popper_index()`](https://nkaza.github.io/gridmorph/reference/gm_polsby_popper_index.md)/etc.'s
convention of omitting the argument entirely rather than silently
ignoring it.

`n_points` points are drawn WITHOUT replacement from the finite set of
boundary cells (capped at however many exist, with a warning if
`n_points` had to be reduced) - unlike
[`gm_geodesic_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_span_index.md)'s
own with-replacement interior sampling, there is no density to weight by
here, so plain uniform sampling over a known finite population is both
simpler and avoids any duplicate-point risk by construction. A
deliberately different argument name from
[`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)'s
own `size` - see
[`gm_geodesic_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_span_index.md)'s
own doc for why.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_geodesic_chord_index(r, n_points = 25, seed = 1)$index
#> [1] 0.9822228
```
