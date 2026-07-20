# Geodesic chord index for a terra raster

Mean GEODESIC distance between two random points on the shape's own
BOUNDARY, relative to an equal-area circle - the raster analogue of
Angel, Parent & Civco (2010)'s "Traversal Index" (see file header for
the naming choice and why it's called "chord" here, not "traversal").
`index = D_ref/D`, in `(0, 1]`, `= 1` iff the shape is (almost
everywhere) a disk.

## Usage

``` r
gm_geodesic_chord_index(rast, size = 40, seed = NULL)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

- size:

  number of boundary points to sample as sources - matches
  [`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)'s
  own argument name and meaning. Checked against a memory/time-derived
  ceiling before running.

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

`size` SOURCE points are drawn WITHOUT replacement from the finite set
of boundary cells (capped at however many exist, with a warning if
`size` had to be reduced) - a genuinely different, and still correct,
convention from
[`gm_geodesic_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_span_index.md)'s
own with-replacement interior sampling: there is no density to weight by
here, so plain uniform sampling over a known finite population is both
simpler and avoids any duplicate-source risk by construction. Each
source's own contribution is the exact UNWEIGHTED mean distance to every
OTHER boundary cell (see
[`gm_geodesic_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_span_index.md)'s
own doc, and the file header, for why this beats sampling a partner).

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_geodesic_chord_index(r, size = 25, seed = 1)$index
#> [1] 0.9823251
```
