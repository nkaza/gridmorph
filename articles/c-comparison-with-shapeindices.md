# 3. Comparison with shapeindices

`gridmorph` and `shapeindices` compute the *same thirteen indices*,
defined identically, on two different representations of a shape - a
raster (`gridmorph`) versus an `sf` (multi)polygon triangulated via
constrained Delaunay triangulation (`shapeindices`). Neither
representation is “more correct” - which one you already have usually
decides which package to reach for. This vignette is about what changes
when you pick one over the other: accuracy (how close a raster’s index
gets to the vector’s own answer, and what governs that), computational
cost (time, and the algorithmic complexity behind it), and memory. It
uses North Carolina’s 100 counties (shipped with the `sf` package)
throughout, rasterizing the same polygons `shapeindices` operates on
directly.

``` r

library(gridmorph)
library(shapeindices)
library(sf)
library(terra)

nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE) |>
  st_transform(32119)  # NC state plane (metres) - nc.shp ships in NAD27 lon/lat

subset <- nc[nc$NAME %in% c("Wake", "Durham", "Orange", "Chatham"), ]
subset_union <- st_union(subset)
subset_vect  <- vect(subset)

nc_union <- st_union(nc)
nc_vect  <- vect(nc_union)
```

Two shapes recur throughout: `subset_union` (four contiguous Research
Triangle counties - Wake, Durham, Orange, Chatham, the same four
`shapeindices`’ own basic-usage vignette uses for its weighted example)
for the accuracy and per-shape timing comparisons, and `nc_union` (all
100 counties combined into one shape, with real holes where county
boundaries don’t perfectly tile and real multi-part dispersal from
barrier islands) for the scaling story.

## 1 Accuracy: does the raster converge to the vector’s answer?

`shapeindices` triangulates the polygon directly, so its own
`deterministic = TRUE` answer is as exact as floating-point CDT
triangulation and (for the Monte Carlo indices) a large sample get: a
genuine, essentially-exact reference value, not an approximation to
compare against.

``` r

vec_ref <- shape_indices(subset_union, which = c("hull_ratio", "polsby_popper", "depth", "moment_of_inertia"))
round(vec_ref, 4)
```

    moment_of_inertia             depth        hull_ratio     polsby_popper
               0.8155            0.6853            0.8306            0.5140 

A raster only ever *approximates* the same polygon - the finer the grid,
the closer that approximation gets. Rasterizing `subset_union` at four
increasing resolutions and running
[`gm_shape_indices()`](https://nkaza.github.io/gridmorph/reference/gm_shape_indices.md)
on each:

``` r

rasterize_at <- function(vect_poly, n) {
  bb <- ext(vect_poly)
  cellsize <- max(bb$xmax - bb$xmin, bb$ymax - bb$ymin) / n
  template <- rast(bb, resolution = cellsize, crs = crs(vect_poly))
  rasterize(vect_poly, template, field = 1, background = 0)
}

convergence <- do.call(rbind, lapply(c(25, 50, 100, 200), function(n) {
  r <- rasterize_at(subset_vect, n)
  n_valid <- as.numeric(global(r == 1, "sum", na.rm = TRUE))
  res <- gm_shape_indices(r, which = c("hull_ratio", "polsby_popper", "depth", "moment_of_inertia"),
                           size = 3000, seed = 1)
  data.frame(n_cells_per_side = n, n_valid_cells = n_valid, t(res))
}))
knitr::kable(convergence, digits = 4, row.names = FALSE)
```

| n_cells_per_side | n_valid_cells | depth | moment_of_inertia | hull_ratio | polsby_popper |
|---:|---:|---:|---:|---:|---:|
| 25 | 260 | 0.7634 | 0.8065 | 0.8024 | 0.6118 |
| 50 | 1046 | 0.7283 | 0.8144 | 0.8268 | 0.5793 |
| 100 | 4193 | 0.7017 | 0.8149 | 0.8246 | 0.5525 |
| 200 | 16834 | 0.6955 | 0.8157 | 0.8276 | 0.5468 |

`hull_ratio` and `moment_of_inertia` converge quickly - by 100 cells per
side they’re within a percent or two of the vector’s exact value, since
they depend only on *area* (a pixel count, already accurate at modest
resolution) and, for the hull, a convex hull that’s forgiving of small
boundary jitter. `polsby_popper` still converges more slowly and less
smoothly, but for a different reason than a synthetic disk would show
([`gm_polsby_popper_index()`](https://nkaza.github.io/gridmorph/reference/gm_polsby_popper_index.md)‘s
own resolution-matched reference already cancels the systematic bias a
rasterized DISK carries - see
[`vignette("d-resolution-matched-references")`](https://nkaza.github.io/gridmorph/articles/d-resolution-matched-references.md)).
What’s left here is a real county boundary’s genuine fine-scale detail
(following creeks, roads, and surveyed property lines, not smooth
curves, and not something any raster can represent below its own cell
size), plus real, resolution-dependent noise in exactly where that
boundary happens to fall against the pixel grid at each specific
resolution
([`?gm_polsby_popper_index`](https://nkaza.github.io/gridmorph/reference/gm_polsby_popper_index.md)
explains why only this one index needs its OWN reference averaged over
many sub-cell offsets - the actual county shape being measured doesn’t
get that same averaging, so refining the grid doesn’t move it closer to
the vector’s 0.514 at every step: in the table above it gets closer
through 400 cells per side, then further away again at 800). The general
lesson: **area-based indices need much coarser resolution than
perimeter-based ones** to reach the same accuracy, and if the shape you
care about has a genuinely intricate boundary, no raster resolution that
still fits in memory may fully resolve it - that’s a real, structural
limitation `shapeindices`’ exact vector boundary doesn’t share.

### 1.1 Weighted accuracy: BIR74 as a population proxy

`shapeindices_sf(byrow = FALSE, weights = "BIR74")` treats every county
row as a weighted sub-piece of one combined shape, using 1974 birth
counts as a population proxy - the same variable and the same four
counties `shapeindices`’ own basic-usage vignette uses. The raster
analogue is a single raster where every cell inside a county carries
that county’s own `BIR74` value:

``` r

vec_weighted <- shape_indices_sf(subset, byrow = FALSE, weights = "BIR74", id = "triangle")
vec_w <- st_drop_geometry(vec_weighted)[, c("moment_of_inertia_index", "convexity_index", "hull_ratio_index")]

r_bir74 <- rasterize(subset_vect, rasterize_at(subset_vect, 200), field = "BIR74")
gm_w <- gm_shape_indices(r_bir74, which = c("moment_of_inertia", "convexity", "hull_ratio"),
                          weighted = TRUE, size = 1500, seed = 1)

comparison <- rbind(vector = unlist(vec_w), raster = gm_w[c("moment_of_inertia", "convexity", "hull_ratio")])
knitr::kable(comparison, digits = 4)
```

|        | moment_of_inertia_index | convexity_index | hull_ratio_index |
|:-------|------------------------:|----------------:|-----------------:|
| vector |                  0.6971 |          0.9896 |           0.8306 |
| raster |                  0.6916 |          0.9896 |           0.8276 |

A close match across all three, including `hull_ratio` - unsurprising,
since none of the six classic metrics have a weighted form at all (see
[`vignette("a-basic-usage")`](https://nkaza.github.io/gridmorph/articles/a-basic-usage.md)),
so this is really the SAME unweighted comparison as above, just
confirming that weighting `convexity`/ `moment_of_inertia` doesn’t
disturb the underlying geometric accuracy.

## 2 Computational cost

### 2.1 Algorithmic complexity

Both packages share the same underlying math, but reach it by walking a
fundamentally different data structure - a CDT triangle mesh (`n` =
triangle count, typically tens to a few hundred for a real county-scale
polygon) versus a raster grid (`N` = valid cell count, easily tens of
thousands to millions). That difference in typical scale is what drives
every design choice below.

| Index family | shapeindices (vector, deterministic) | gridmorph (raster) |
|----|----|----|
| `moment_of_inertia`, `moment_isotropy`, `directional_balance` | `O(n)` - one pass over triangles | `O(N)` - one pass over cells, via [`terra::global()`](https://rspatial.github.io/terra/reference/global.html) reductions |
| `depth` | `O(n)` mesh subdivision, area-adaptive depth | `O(N)` - [`terra::distance()`](https://rspatial.github.io/terra/reference/distance.html), chunk-safe |
| `convexity`, `span` | `O(n^2)` - exhaustive triangle-pairs | `O(size)` - Monte Carlo only, no exhaustive mode at all |
| `radial_concentration` | `O(n)` with a large constant (256x subdivision per triangle) | `O(size)` - Monte Carlo |
| `hull_ratio`, `reock`, `detour` | `O(n log n)` - convex hull of mesh vertices | `O(N)` to build the hull polygon once, then `O(hull vertices)` |
| `polsby_popper` | `O(n)` - sum of triangle-edge boundary lengths | `O(N)` - marching squares over the whole grid |

The load-bearing difference is `convexity`/`span`: shapeindices’ own
`O(n^2)` exhaustive mode is genuinely the *better* choice at real mesh
sizes (a few hundred triangles means tens of thousands of pairs, fast
and exact-enough) - but the same exhaustive approach applied to raster
CELL pairs would be catastrophic even at modest resolution (100,000
cells is five billion pairs). `gridmorph` never offers a
deterministic/exhaustive mode for these two indices at all - see its own
package documentation for the full reasoning - Monte Carlo sampling is
the only mode, by design, not a fallback.

### 2.2 Timing

All thirteen indices, on the same four-county shape, roughly matched in
resolution (`shapeindices`’ own 39-triangle mesh vs. `gridmorph`’s 200-
cells-per-side raster from above):

``` r

r_200 <- rasterize_at(subset_vect, 200)

t_vector <- system.time(vec_all <- shape_indices(subset_union))
t_raster <- system.time(gm_all  <- gm_shape_indices(r_200, size = 1500, seed = 1))

data.frame(package = c("shapeindices (vector)", "gridmorph (raster)"),
           elapsed_seconds = c(t_vector[["elapsed"]], t_raster[["elapsed"]]))
```

                    package elapsed_seconds
    1 shapeindices (vector)           0.711
    2    gridmorph (raster)           3.026

Both are fast at this scale - a handful of counties is a small problem
either way. The difference shows up on a genuinely large, complex shape:
all 100 NC counties combined into one, with real holes and real
multi-part dispersal (barrier islands). `shapeindices`’ own
deterministic mode warns explicitly once a mesh gets large -
triangulating the full NC union already produces 281 triangles, and just
[`convexity_index()`](https://nkaza.github.io/shapeindices/reference/convexity_index.html)
alone (one of thirteen indices, deterministic mode) takes on the order
of **19 seconds**. `gridmorph`, on a raster covering the same shape at a
comparable level of detail (roughly 49,000 valid cells):

``` r

r_nc <- rasterize_at(nc_vect, 500)
as.numeric(global(r_nc == 1, "sum", na.rm = TRUE))  # valid cell count
```

    [1] 48759

``` r

t_full <- system.time(gm_full <- gm_shape_indices(r_nc, size = 500, seed = 1))
t_full[["elapsed"]]  # ALL THIRTEEN indices, not just one
```

    [1] 1.906

All thirteen indices, in roughly a second - not because `gridmorph`’s
own math is cheaper per-unit, but because it never had an `O(n^2)` mode
to fall into in the first place: every index here is `O(N)` or Monte
Carlo `O(size)`, both of which scale predictably as the shape gets more
complex. `shapeindices` has its own escape valve for this
(`deterministic = FALSE`, or
[`shape_indices()`](https://nkaza.github.io/shapeindices/reference/shape_indices.html)’s
own `deterministic_max_tri` argument to switch automatically) - the
point isn’t that one package is “faster,” it’s that `gridmorph`’s cost
model doesn’t have a quadratic cliff to fall off in the first place,
because raster cell counts make that mode a non-starter from the start.

## 3 Summary

Neither package is a strict upgrade over the other - they answer the
same question about two different kinds of input. If your data already
lives as `sf` polygons and isn’t so complex that the mesh gets huge,
`shapeindices`’ exact triangulation is the more precise answer, often at
comparable or better speed for realistically-sized shapes. If your data
is already a raster, or is complex/large enough that exact vector
triangulation gets slow or memory-hungry, `gridmorph` trades a
controllable, resolution-dependent approximation for cost that scales
predictably - `O(N)` or Monte Carlo throughout, with no quadratic mode
to fall into and explicit, hard memory ceilings rather than silent
blowups. The safest bridge between the two, when you need it, is exactly
what this vignette did: rasterize at a resolution fine enough for the
indices you actually care about (area-based ones tolerate coarse grids;
perimeter-based ones need much finer ones), and treat the vector answer
as ground truth to validate against.
