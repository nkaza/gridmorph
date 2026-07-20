# Morphological dilation of a terra SpatRaster

For binary (0/1) or continuous input, the local maximum over `kernel`'s
footprint at every cell (grayscale dilation, a strict generalization of
binary dilation, not a different operation); for a categorical
([`terra::is.factor()`](https://rspatial.github.io/terra/reference/is.bool.html))
input, standard label-image dilation - already-labelled cells are left
untouched, and `NA` cells are filled with the majority label among their
real-valued neighbours (see this file's own header for the full
reasoning behind both).

## Usage

``` r
dilate(mask, kernel = se_box(1))
```

## Arguments

- mask:

  binary, continuous, or categorical SpatRaster

- kernel:

  a structuring element: any matrix, or one of
  [`se_box()`](https://nkaza.github.io/gridmorph/reference/se_box.md),
  [`se_disc()`](https://nkaza.github.io/gridmorph/reference/se_disc.md),
  [`se_diamond()`](https://nkaza.github.io/gridmorph/reference/se_diamond.md).
  `1` includes a cell in the footprint, `NA` excludes it (matching
  [`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)'s
  own `w` convention) - any other value errors (see `.check_kernel()`'s
  own comments for why).

## Value

SpatRaster, same grid and value type as `mask` (see
[`erode()`](https://nkaza.github.io/gridmorph/reference/erode.md))

## Examples

``` r
r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, crs = "local")
terra::values(r) <- 0
r[8:13, 8:13] <- 1
dilate(r, kernel = se_disc(1))
#> class       : SpatRaster
#> size        : 20, 20, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 20, 0, 20  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter)
#> source(s)   : memory
#> name        : focal_max
#> min value   :         0
#> max value   :         1
```
