# Top-hat (white-hat) transform of a terra SpatRaster

Highlights small bright features/protrusions narrower than the
structuring element. For binary/continuous input, the real-valued
residual `mask - opening(mask)` (reduces exactly to the classic
`mask AND NOT opening(mask)` when input is 0/1). For categorical input,
a boolean indicator instead - `1` at cells that had a label but lost it
under
[`opening()`](https://nkaza.github.io/gridmorph/reference/opening.md) -
since subtracting category codes is meaningless (see this file's own
header for why).

## Usage

``` r
tophat(mask, kernel = se_box(1))
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

SpatRaster, same grid as `mask`. Binary/continuous input: `NA` outside
`mask`'s own footprint, non-negative and real-valued inside it wherever
`kernel` includes its own centre (true for
[`se_box()`](https://nkaza.github.io/gridmorph/reference/se_box.md)/[`se_disc()`](https://nkaza.github.io/gridmorph/reference/se_disc.md)/[`se_diamond()`](https://nkaza.github.io/gridmorph/reference/se_diamond.md)
at every radius - a custom kernel that excludes its own centre can
violate this, a property of grayscale morphology in general).
Categorical input: plain `0`/`1`, not itself categorical - the result is
a flag, not a category.

## Details

Binary/continuous input specifically: `NA` outside `mask`'s own
footprint, not a numeric `0` -
[`opening()`](https://nkaza.github.io/gridmorph/reference/opening.md)'s
own output never actually contains `NA` for non-categorical input
(`.mask01()` replaces `NA` with `0` internally so `focal()` treats it as
background), so without this restriction the residual would be a
fully-defined `0` arbitrarily far outside the shape. Restoring `NA`
there loses no information: `mask - opening(mask)` is exactly `0`
everywhere outside `mask`'s own footprint regardless (`opening(x) <= x`
pointwise, and opening never creates value beyond `mask`'s own extent).
