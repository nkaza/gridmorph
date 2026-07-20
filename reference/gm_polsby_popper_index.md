# Polsby-Popper compactness score of a terra raster

`raw / ref_index`, in `(0, 1]`, `1` = a circle, where `raw` is the
textbook `4*pi*area/perimeter^2` (perimeter via marching squares,
[`terra::as.contour()`](https://rspatial.github.io/terra/reference/contour.html) -
see this file's own header for why, verified empirically against the
naive cell-edge alternative) and `ref_index` is that SAME formula's own
score on an actual rasterized disk at `rast`'s own cell size and area.

## Usage

``` r
gm_polsby_popper_index(rast)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

## Value

`list(index, area, perimeter, ref_index, n_valid_cells)`. `ref_index` is
the reference disk's own raw score - `index` itself is
`(4*pi*area/perimeter^2) / ref_index`, not the raw formula directly.

## Details

WHY NOT JUST `raw`: comparing it against the textbook constant `1` (a
claim about a perfect, infinitely-smooth circle) systematically
understates how compact a rasterized shape is. A rasterized DISK - the
one shape this score should never meaningfully penalize - scored `0.88`
on `raw` alone at a representative resolution, not because the shape is
imperfect, but because marching-squares perimeter overestimates a raster
boundary's true length, and squaring it in the denominator roughly
doubles that percentage error. Dividing by `ref_index` instead cancels
that measurement bias, since both numerator and reference go through the
identical perimeter-measurement pipeline: the same disk scores
`0.97-1.00` this way, across resolutions. This reference is
PHASE-AVERAGED (25 sub-cell offsets) rather than a single realization -
see this file's own header for why marching-squares perimeter needs that
where the other resolution-matched references in this file don't. No
`analytical_ref` opt-out is provided: the raw closed-form score is not a
more "correct" number to fall back to, it is measurably the wrong one
for a rasterized shape at any finite resolution. No `weighted`
argument - see this file's own header for why.

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_polsby_popper_index(r)$index
#> [1] 0.8989195
```
