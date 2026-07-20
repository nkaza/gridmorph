# All indices (or a chosen subset) for a single terra raster

All indices (or a chosen subset) for a single terra raster

## Usage

``` r
gm_shape_indices(rast, which = "all", ...)
```

## Arguments

- rast:

  a terra SpatRaster - see
  [`gm_depth_index()`](https://nkaza.github.io/gridmorph/reference/gm_depth_index.md)
  for the shared shape/hole conventions every index in this package
  uses.

- which:

  `"all"` (default, expands to all fifteen values below), or a character
  vector naming a subset - each listed here with the function it
  actually calls, since the `which` string and the function name aren't
  identical:

  - `"depth"` -
    [`gm_depth_index()`](https://nkaza.github.io/gridmorph/reference/gm_depth_index.md)

  - `"moment_of_inertia"` -
    [`gm_moment_of_inertia_index()`](https://nkaza.github.io/gridmorph/reference/gm_moment_of_inertia_index.md)

  - `"moment_isotropy"` -
    [`gm_moment_isotropy_index()`](https://nkaza.github.io/gridmorph/reference/gm_moment_isotropy_index.md)

  - `"directional_balance"` -
    [`gm_directional_balance_index()`](https://nkaza.github.io/gridmorph/reference/gm_directional_balance_index.md)

  - `"convexity"` -
    [`gm_convexity_index()`](https://nkaza.github.io/gridmorph/reference/gm_convexity_index.md)

  - `"span"` -
    [`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)

  - `"radial_concentration"` -
    [`gm_radial_concentration_index()`](https://nkaza.github.io/gridmorph/reference/gm_radial_concentration_index.md)

  - `"hull_ratio"` -
    [`gm_hull_ratio_index()`](https://nkaza.github.io/gridmorph/reference/gm_hull_ratio_index.md)

  - `"polsby_popper"` -
    [`gm_polsby_popper_index()`](https://nkaza.github.io/gridmorph/reference/gm_polsby_popper_index.md)

  - `"width_length_ratio"` -
    [`gm_width_length_ratio_index()`](https://nkaza.github.io/gridmorph/reference/gm_width_length_ratio_index.md)

  - `"reock"` -
    [`gm_reock_index()`](https://nkaza.github.io/gridmorph/reference/gm_reock_index.md)

  - `"detour"` -
    [`gm_detour_index()`](https://nkaza.github.io/gridmorph/reference/gm_detour_index.md)

  - `"exchange"` -
    [`gm_exchange_index()`](https://nkaza.github.io/gridmorph/reference/gm_exchange_index.md)

  - `"geodesic_span"` -
    [`gm_geodesic_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_span_index.md) -
    substantially more expensive than every index above (`n_points`
    sequential whole-raster
    [`terra::gridDist()`](https://rspatial.github.io/terra/reference/gridDist.html)
    calls, not `O(size)` or cheaper - see that function's own file
    header) but included in `"all"` regardless, on the reasoning that
    `"all"` should mean all fifteen, not
    thirteen-plus-two-you-have-to-know-to-ask-for. Its own sample size
    is `n_points`, a DELIBERATELY DIFFERENT argument from `size` (see
    its own doc) - lower `n_points` (via `...`) for a faster `"all"`
    call if this matters; passing `size` here has no effect on it at
    all.

  - `"geodesic_chord"` -
    [`gm_geodesic_chord_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_chord_index.md) -
    same cost note and same `n_points` (not `size`) argument as
    `"geodesic_span"` above.

  The first thirteen names are the same short names
  [`shapeindices::shape_indices()`](https://nkaza.github.io/shapeindices/reference/shape_indices.html)
  uses, so results from the two packages line up directly for those
  (`shapeindices` has no geodesic-distance indices to compare against -
  see `shapeindices/explorations/NOTES-future-indices.md`). An
  unrecognised name in `which` errors immediately, listing all fifteen
  valid values.

- ...:

  passed to whichever of the requested index functions accept each named
  argument (e.g. `weighted` is accepted by `"depth"` through
  `"radial_concentration"` and `"geodesic_span"`, silently ignored for
  the six classic metrics and `"geodesic_chord"`, none of which have a
  weighted form - see each's own file header for why; `size`/`seed` are
  accepted by the three ORIGINAL Monte Carlo indices
  ([`gm_convexity_index()`](https://nkaza.github.io/gridmorph/reference/gm_convexity_index.md)/[`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)/[`gm_radial_concentration_index()`](https://nkaza.github.io/gridmorph/reference/gm_radial_concentration_index.md));
  `n_points`/`seed` are accepted by
  [`gm_geodesic_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_span_index.md)/
  [`gm_geodesic_chord_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_chord_index.md)
  instead - a DELIBERATELY DIFFERENT argument name from `size`, not an
  inconsistency: these two cost `O(n_points * n_cells)`, not `O(size)`,
  so sharing one argument name would mean an ordinary `size = 3000` call
  (sensible for the other three) silently driving 3000 sequential
  whole-raster
  [`terra::gridDist()`](https://rspatial.github.io/terra/reference/gridDist.html)
  calls too - verified directly to cause a real, surprising slowdown
  before this split (see R/geodesic-index.R's own file header); `n_bins`
  is accepted by
  [`gm_depth_index()`](https://nkaza.github.io/gridmorph/reference/gm_depth_index.md)/
  [`gm_moment_of_inertia_index()`](https://nkaza.github.io/gridmorph/reference/gm_moment_of_inertia_index.md)/[`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)/
  [`gm_radial_concentration_index()`](https://nkaza.github.io/gridmorph/reference/gm_radial_concentration_index.md)/[`gm_geodesic_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_geodesic_span_index.md) -
  passing it explicitly overrides each of their own individual defaults,
  including
  [`gm_span_index()`](https://nkaza.github.io/gridmorph/reference/gm_span_index.md)'s
  deliberately smaller one).

## Value

named numeric vector, one entry per requested index, in canonical order

## Examples

``` r
r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "local")
terra::values(r) <- 0
r[10:30, 10:30] <- 1
gm_shape_indices(r, size = 1000, seed = 1)
#>                depth    moment_of_inertia      moment_isotropy 
#>            0.9340579            0.9571000            1.0000000 
#>  directional_balance            convexity                 span 
#>            0.9977324            1.0000000            0.9754981 
#> radial_concentration           hull_ratio        polsby_popper 
#>            0.9957273            1.0000000            0.8989195 
#>   width_length_ratio                reock               detour 
#>            1.0000000            0.6979405            0.9150637 
#>             exchange        geodesic_span       geodesic_chord 
#>            0.9092971            0.8488987            0.9686795 

# a subset
gm_shape_indices(r, which = c("hull_ratio", "polsby_popper", "reock"))
#>    hull_ratio polsby_popper         reock 
#>     1.0000000     0.8989195     0.6979405 
```
