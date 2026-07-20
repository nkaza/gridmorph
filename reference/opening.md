# Morphological opening (erode then dilate) of a terra SpatRaster

Removes small bright features and thin protrusions without shrinking the
overall shape the way
[`erode()`](https://nkaza.github.io/gridmorph/reference/erode.md) alone
does. Generalizes to continuous input as grayscale opening, and to
categorical input as label-image opening (removes small/thin
single-class regions), the same operation in each case, via
[`erode()`](https://nkaza.github.io/gridmorph/reference/erode.md)/[`dilate()`](https://nkaza.github.io/gridmorph/reference/dilate.md)'s
own dispatch.

## Usage

``` r
opening(mask, kernel = se_box(1))
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
