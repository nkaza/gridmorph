# Moment isotropy index for a terra raster

The raster analogue of
[`shapeindices::moment_isotropy_index()`](https://nkaza.github.io/shapeindices/reference/moment_isotropy_index.html).
index in `(0, 1]`, ratio of the smaller to larger principal moment of
the mass inertia tensor - `= 1` iff the mass distribution is
rotationally isotropic about its own centroid (any shape with 3-fold or
higher rotational symmetry qualifies, not only a disk).

## Usage

``` r
gm_moment_isotropy_index(rast, weighted = TRUE)
```

## Arguments

- rast:

  a terra SpatRaster. The shape is derived directly from `rast`: a cell
  is part of the shape iff its own value is neither `NA` nor exactly
  `0` - both are holes, no separate mask argument.

- weighted:

  logical, default `TRUE`. When `TRUE`, `rast`'s own cell values are
  used as the density/mass throughout. When `FALSE`, every valid cell is
  treated as equally massed regardless of `rast`'s actual values -
  exactly reproducing the plain (unweighted) index even if `rast` is
  itself a continuous raster. If `rast` is a CATEGORICAL raster
  (`terra::is.factor(rast)` is `TRUE` - e.g. a land-use classification
  with an attached levels table) and `weighted = TRUE`, the category
  CODE NUMBERS themselves get used as mass, which is rarely the intended
  comparison (a class coded `3` outweighs one coded `1` for no
  meaningful reason) - a warning is issued but the computation still
  proceeds on the literal codes, since `weighted = TRUE` was explicitly
  requested; pass `weighted = FALSE` for categorical rasters unless the
  codes really are meant as magnitudes. This check never fires on
  non-categorical (including plain binary 0/1) rasters, so
  `weighted = TRUE` is always safe there.

## Value

list(index, Ixx, Iyy, Ixy, centroid)
