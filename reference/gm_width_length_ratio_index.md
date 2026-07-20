# Width-length ratio of a terra raster's bounding box

The shorter of the bounding box's x/y extents over the longer, in
`(0, 1]`, `1` = a square bounding box. Axis-aligned, not the minimum
bounding rectangle at any rotation - a diagonally-oriented elongated
shape can score deceptively high, the classic limitation of this score.
No `weighted` argument - see this file's own header for why.

## Usage

``` r
gm_width_length_ratio_index(rast)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

## Value

`list(index, length, width, n_valid_cells)`

## Details

KNOWN LIMITATION (ported unchanged from shapeindices): blind to both
holes and multi-part dispersal, since only the bounding box's own extent
enters the ratio - a shape and the same shape with a large hole punched
through it score identically.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_width_length_ratio_index(r)$index
#> [1] 1
```
