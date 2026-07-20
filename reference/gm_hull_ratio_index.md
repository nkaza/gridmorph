# Convex hull area ratio of a terra raster

area / area(convex hull), in `(0, 1]`, `1` = the shape is convex (equals
its own convex hull). No `weighted` argument - see this file's own
header for why area/hull are structurally unweighted quantities.

## Usage

``` r
gm_hull_ratio_index(rast)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

## Value

`list(index, area, hull_area, hull, n_valid_cells)`. `hull` is the
convex hull itself, a terra SpatVector, for plotting.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_hull_ratio_index(r)$index
#> [1] 1
```
