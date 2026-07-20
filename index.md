Computes shape/compactness/dispersal indices for a single terra
SpatRaster, treating a whole shape (however many disjoint patches or
holes it has) as a single combined shape - the raster- native sibling of
'shapeindices', which does the same thing for (multi)polygons via
constrained Delaunay triangulation. The shape is derived directly from
the raster's own values, optionally used as a density/mass field
throughout. Every index requires a planar (projected) coordinate
reference system, since the package's own point-based geometry
(line-crossing tests, pairwise distance accumulation, moment-tensor
calculations, circle-fitting) works directly in the raster's own x/y
coordinate space as if it were Cartesian - a geographic
(longitude/latitude) CRS is rejected with an error rather than silently
producing wrong results. Also exposes morphological operators (erosion,
dilation, opening, closing, top-hat, bottom-hat) as first-class
functions on 'terra' SpatRaster objects, filling a gap 'terra' itself
does not cover.
