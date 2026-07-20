# Reock compactness score of a terra raster

`raw / ref_index`, in `(0, 1]`, `1` = a circle, where `raw` is
`area/area(minimum bounding circle)` and `ref_index` is that SAME
formula's own score on an actual rasterized disk at `rast`'s own cell
size and area. The minimum bounding circle depends only on the shape's
own convex hull vertices, found via a deterministic (no RNG) Welzl-style
algorithm - exact, not an approximation.

## Usage

``` r
gm_reock_index(rast)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

## Value

`list(index, area, mbc_area, ref_index, mbc, n_valid_cells)`. `mbc` is
the minimum bounding circle itself, a terra SpatVector, for plotting.
`ref_index` is the reference disk's own raw score - `index` itself is
`(area/mbc_area) / ref_index`, not the raw ratio directly.

## Details

WHY NOT JUST `raw`: `area` is a pixel COUNT, which undershoots a shape's
true continuous area near its own boundary (a boundary cell only counts
as "in" when its centre does), while `mbc_area` comes from EXACT hull
vertices - two different measurement conventions for what should be the
same disk. A rasterized DISK scored `raw = 0.953` at a representative
resolution, not because the shape is imperfect, but from that convention
mismatch alone. Dividing by `ref_index` instead cancels it, since both
numerator and reference go through the identical pixel-count-area /
exact-MBC-area pipeline: the same disk scores `0.996-1.000` this way,
across resolutions (verified: low phase sensitivity here, unlike
[`gm_polsby_popper_index()`](https://nkaza.github.io/gridmorph/reference/gm_polsby_popper_index.md)'s
own marching- squares perimeter - see this file's own header). No
`analytical_ref` opt-out is provided: the raw ratio is not a more
"correct" number to fall back to, it is measurably the wrong one for a
rasterized shape at any finite resolution. No `weighted` argument - see
this file's own header for why.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_reock_index(r)$index
#> [1] 0.6979405
```
