# Width-length ratio of a terra raster's minimum-area bounding rectangle

The shorter side of the shape's minimum-area bounding rectangle, at ANY
rotation, over the longer side, in `(0, 1]`, `1` = that rectangle is a
square. Found via rotating calipers over the shape's own convex hull
([`terra::convHull()`](https://rspatial.github.io/terra/reference/convhull.html)) -
deliberately NOT the raster's axis-aligned extent, matching a fix
already made to
[`shapeindices::width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.html):
an axis-aligned box scores a shape's elongation relative to how it
happens to be oriented on the grid, not relative to the shape itself -
rotating a shape in place, without changing it at all otherwise, could
swing the old score anywhere from its true value up to a coincidental
`1`. No `weighted` argument - see this file's own header for why.

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

`list(index, length, width, mbr, n_valid_cells)`. `mbr` is the
minimum-area bounding rectangle itself, a terra SpatVector, for
plotting.

## Details

KNOWN LIMITATION (ported unchanged from shapeindices): blind to both
holes and multi-part dispersal, since only the bounding rectangle's own
extent enters the ratio - a shape and the same shape with a large hole
punched through it score identically.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_width_length_ratio_index(r)$index
#> [1] 1
```
