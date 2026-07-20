# Detour compactness score of a terra raster

`raw / ref_index`, in `(0, 1]`, `1` = a circle, where `raw` is the ratio
of the equal-area circle's ANALYTICAL perimeter to the perimeter of the
shape's own convex hull, and `ref_index` is that SAME formula's own
score on an actual rasterized disk at `rast`'s own cell size and area.
Angel, Parent & Civco (2010) introduce this to measure how hard a shape
is to circumnavigate as an obstacle. The hull's own perimeter needs no
marching-squares correction (unlike
[`gm_polsby_popper_index()`](https://nkaza.github.io/gridmorph/reference/gm_polsby_popper_index.md)'s
raw shape perimeter) -
[`terra::convHull()`](https://rspatial.github.io/terra/reference/convhull.html)'s
output is already a proper straight-edged polygon, not a raster-traced
boundary.

## Usage

``` r
gm_detour_index(rast)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

## Value

`list(index, area, hull_perimeter, ref_index, hull, n_valid_cells)`.
`hull` is the convex hull itself, a terra SpatVector, for plotting.
`ref_index` is the reference disk's own raw score - `index` itself is
`(2*sqrt(pi*area)/hull_perimeter) / ref_index`, not the raw ratio
directly.

## Details

WHY NOT JUST `raw`: under centre-based rasterization, a boundary cell's
own far CORNER can sit outside the shape's true continuous boundary by
up to half a cell's diagonal, so the hull built from those corners has
vertices that genuinely protrude past the true circle - its perimeter
comes out LARGER than the analytical circle's own, even for an actual
disk. That scored `raw = 0.982` at a representative resolution. Dividing
by `ref_index` instead cancels it, since both numerator and reference go
through the identical hull-perimeter pipeline: the same disk scores
`0.997-1.000` this way, across resolutions (verified: low phase
sensitivity here, unlike
[`gm_polsby_popper_index()`](https://nkaza.github.io/gridmorph/reference/gm_polsby_popper_index.md)'s
own marching-squares perimeter - see this file's own header). No
`analytical_ref` opt-out is provided: the raw ratio is not a more
"correct" number to fall back to, it is measurably the wrong one for a
rasterized shape at any finite resolution. No `weighted` argument - see
this file's own header for why.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_detour_index(r)$index
#> [1] 0.9150637
```
