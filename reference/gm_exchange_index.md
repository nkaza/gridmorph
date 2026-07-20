# Exchange compactness score of a terra raster

The share of the shape's own area that falls inside the equal-area
circle centred at its (unweighted) centroid, in `[0, 1]`, `1` = a
circle. Angel, Parent & Civco (2010) introduce this as a natural metric
for gerrymandering. Computed via pixel counting - each valid cell's own
centre is tested against the reference circle directly, the
raster-native analogue of the vector package's exact polygon
intersection. No `weighted` argument - see this file's own header for
why.

## Usage

``` r
gm_exchange_index(rast)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

## Value

`list(index, area, circle_area, circle, n_valid_cells)`. `circle` is the
equal-area reference circle itself, a terra SpatVector, for plotting.

## Details

KNOWN LIMITATION (ported unchanged from shapeindices): for a multi-part
shape, the reference circle is centred at the OVERALL centroid, which
for well-separated parts can sit in the empty space between them - once
the circle's radius is smaller than the distance to the nearest part,
the index is exactly 0, not just low.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_exchange_index(r)$index
#> [1] 0.9092971
```
