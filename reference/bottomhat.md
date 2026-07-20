# Bottom-hat (black-hat) transform of a terra SpatRaster

Highlights small dark gaps/notches narrower than the structuring
element. For binary/continuous input, the real-valued residual
`closing(mask) - mask` (reduces exactly to the classic
`closing(mask) AND NOT mask` when input is 0/1). For categorical input,
a boolean indicator instead - `1` at cells that were unlabelled (`NA`)
but got filled by
[`closing()`](https://nkaza.github.io/gridmorph/reference/closing.md) -
since subtracting category codes is meaningless (see this file's own
header for why).

## Usage

``` r
bottomhat(mask, kernel = se_box(1))
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

SpatRaster, same grid as `mask` (see
[`tophat()`](https://nkaza.github.io/gridmorph/reference/tophat.md)'s
own note on value type and non-negativity; `NA` region described above,
not the same one)

## Details

Binary/continuous input specifically: `NA` only where `mask` was ITSELF
`NA` (genuinely missing data, not just background) AND beyond one
`kernel` radius of `mask`'s own footprint - narrower than
[`tophat()`](https://nkaza.github.io/gridmorph/reference/tophat.md)'s
equivalent restriction, deliberately, on both counts. Not "wherever
`mask` is `NA`" the way
[`tophat()`](https://nkaza.github.io/gridmorph/reference/tophat.md) is:
`bottomhat()`'s entire purpose is flagging small gaps ADJACENT to (not
inside) the shape, and a gap IS an `NA` cell in `mask` - restoring `NA`
there would erase exactly the cells this function exists to highlight.
And not "wherever `mask` is confirmed `0`" either: a plain binary
(`0`/`1`, no `NA` at all) mask has no genuinely missing data anywhere,
so its `bottomhat()` output stays `0`/`1` everywhere, matching the
classic definition exactly. The one-kernel-radius reach bound is exact,
not a heuristic margin:
[`closing()`](https://nkaza.github.io/gridmorph/reference/closing.md)
can only ever differ from background within one kernel radius of
`mask`'s own footprint (dilate grows by at most one radius, the
subsequent erode can only give at most that much back).
