# Morphological closing (dilate then erode) of a terra SpatRaster

Fills small dark gaps/notches without growing the overall shape the way
[`dilate()`](https://nkaza.github.io/gridmorph/reference/dilate.md)
alone does. Generalizes to continuous input as grayscale closing, and to
categorical input as label-image closing (fills small gaps within a
class region), the same operation in each case, via
[`erode()`](https://nkaza.github.io/gridmorph/reference/erode.md)/[`dilate()`](https://nkaza.github.io/gridmorph/reference/dilate.md)'s
own dispatch.

## Usage

``` r
closing(mask, kernel = se_box(1))
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
