## =========================================================================
## Morphological operators on terra SpatRaster masks - erode, dilate,
## opening, closing, tophat, bottomhat. First-class exports, not an
## internal helper: `terra` has no erosion/dilation with configurable
## structuring elements at all, and `gm_depth_index()` (elsewhere in this
## package) is built on the same "distance field over a mask" idea.
##
## BUILT ON terra::focal(), NOT mmand - a real, verified architecture
## decision, not the original one. mmand::erode()/dilate() are faster
## per-call at sizes that already fit in memory (verified: ~25% faster for
## a box kernel, ~3x faster for a disc kernel, both on a 1M-cell raster),
## but mmand operates on an in-memory array via as.matrix(), full stop -
## no way to run it on a raster too large to fit in memory. terra::focal()
## is chunk-aware on disk-backed SpatRasters the same way
## terra::distance()/patches() already are, which matters more than the
## per-call speed difference for a package meant to handle rasters that
## may not fit in memory. Verified terra::focal(fun=min/max) reproduces
## mmand's own erode()/dilate() output exactly (identical() check, both
## box and disc kernels) before making this switch, not just assumed.
##
## STRUCTURING ELEMENTS: the `kernel` argument to every function here is a
## plain matrix - any matrix the caller builds, or one of the named
## shortcuts below (`se_box()`, `se_disc()`, `se_diamond()`). This matches
## terra::focal()'s own `w` argument directly (no conversion layer needed):
## a `1` includes a cell in the structuring element's footprint; `NA`
## excludes it (verified: focal(fun=min, na.rm=TRUE) with an NA-masked
## weight matrix correctly reproduces mmand's own non-rectangular kernel
## behaviour, e.g. a disc, not just a box) - this holds regardless of
## which `fun` is used (min/max/modal), verified directly for all three.
##
## KERNEL VALIDATION: every function here rejects a kernel whose non-NA
## values aren't all exactly `1` - not a style preference, a real
## correctness gap closed after verifying it directly. `terra::focal()`
## treats a non-0/1/NA `w` matrix as MULTIPLICATIVE weights, not a
## boolean include/exclude mask: a kernel entry of `0.5` made `erode()`
## return `0.5` at a cell that should have stayed `1`; a kernel entry of
## `2` made `dilate()` return `2`; and using `0` (instead of `NA`) to mean
## "exclude this position" made `erode()` return an all-zero raster with
## NO error, because `0 * anything = 0` always wins the `min`.
## `.check_kernel()` catches all of these up front with a clear error
## instead of a silently corrupted result.
##
## RASTER-EDGE PADDING - two separate, verified terra quirks, not the same
## issue, both concerning what happens where `kernel`'s footprint extends
## past the raster's own true extent:
##
## 1. `terra::focal()` always pads beyond the raster's edge (there is no
##    way to turn padding off, only to choose its value via `fillvalue =`
##    or use `expand = TRUE` to replicate the nearest real cell instead -
##    see rspatial/terra#243, a real upstream issue about exactly this,
##    where the default combination of NA padding + na.rm = TRUE silently
##    produced wrong Sobel-filter output at raster edges for the same
##    underlying reason). For MIN (erosion), that default combination is a
##    real bug here too: na.rm = TRUE DROPS the padded (NA) neighbours
##    from the calculation instead of counting them as background, so a
##    shape reaching the raster's own edge failed to erode there at all -
##    verified directly (a 5x5 all-foreground raster eroded by a 3x3 box
##    came back completely unchanged; the correct answer shrinks it to the
##    interior 3x3 block). Fixed in `.erode_core()` with `fillvalue = 0`,
##    consistent with `.mask01()`'s own "NA/edge = background = 0"
##    convention for interior holes. For MAX (dilation), the SAME default
##    combination is NOT a bug - `0` is `max()`'s identity element, so
##    dropping the padding or including it as `0` give an identical
##    result - verified directly, not just assumed symmetric with erode.
##
## 2. `terra::ifel(cond, yes, no)`, before terra 1.9-11, leaked a raster
##    branch's own value at an NA cell in `cond` instead of giving `NA`,
##    whenever `yes` or `no` was itself a SpatRaster rather than a plain
##    scalar (rspatial/terra#2058; fixed upstream, see terra's own
##    NEWS.md). DESCRIPTION requires `terra (>= 1.9-11)` for exactly this
##    reason. `.erode_core_categorical()`'s condition is still built from
##    `is.na()` rather than `==` as cheap defense-in-depth regardless: `is.na(lo)
##    | is.na(hi) | (lo != hi)` can never itself be `NA` by construction,
##    so there's no NA condition cell left for `ifel()` to mishandle.
##
## INPUT VALUE TYPE - three genuinely different cases, dispatched on
## `terra::is.factor(mask)` in `.erode_core()`/`.dilate_core()` (so
## `opening()`/`closing()`/`tophat()`/`bottomhat()`, all built from those
## two, get the right behaviour automatically at every stage - including
## an intermediate result produced by one of them, since the dispatch
## re-checks `is.factor()` on whatever raster it's actually given, not
## just on the original top-level input):
##
## 1. BINARY (0/1, or 0-and-NA/nonzero) or CONTINUOUS, non-factor, input:
##    `erode()`/`dilate()` are `focal(fun = min/max)` with no binarization
##    step (`.mask01()` only ever replaces `NA` with `0`) - this IS the
##    textbook definition of GRAYSCALE morphological erosion/dilation (a
##    local min/max filter), a strict generalization of binary erosion/
##    dilation, not a different operation. `opening()`/`closing()`
##    generalize for the same reason. `tophat()`/`bottomhat()` are the
##    real grayscale residuals - `mask - opening(mask)` and
##    `closing(mask) - mask` - which reduce EXACTLY (not approximately) to
##    the classic binary `mask AND NOT opening(mask)` /
##    `closing(mask) AND NOT mask` when input is 0/1, since
##    `opening(x) <= x` and `closing(x) >= x` pointwise whenever the
##    structuring element includes its own centre (true for every one of
##    `se_box()`/`se_disc()`/`se_diamond()`, at every `radius >= 0`).
##
## 2. CATEGORICAL (`terra::is.factor(mask)`) input - genuinely different
##    operations, not the same min/max/subtraction formulas applied to
##    numeric category codes (which would be arithmetic over an arbitrary
##    coding scheme, meaningless for nominal classes - the exact failure
##    mode this design avoids). Standard label-image morphology instead:
##
##    - `erode()`: a cell KEEPS its own label only if every cell in its
##      neighbourhood shares that EXACT label; otherwise it becomes `NA`.
##      Equivalent to eroding each class's own binary mask separately and
##      taking the union - classes are mutually exclusive, so there's no
##      overlap to resolve - implemented as one pass, not a per-class
##      loop: `focal(fun = min) == focal(fun = max)` over the raw integer
##      codes is true at a cell iff every value in its neighbourhood is
##      identical (for real numbers, min == max implies uniform), and
##      since the neighbourhood always includes the centre itself, that
##      shared value must be the centre's own label. `na.rm = FALSE`
##      deliberately (not `TRUE`, unlike every other focal() call in this
##      file): an `NA` neighbour must break unanimity and erode the cell,
##      the same "a hole erodes what's next to it" behaviour `.mask01()`
##      already gives the binary/continuous case via forcing NA to `0`
##      (the local minimum) - here there's no single "background level"
##      to fall back on, so the NA has to propagate through the
##      min/max comparison directly instead. `terra::focal()`'s numeric
##      min/max strip the `is.factor()`/`levels()` attribute even when
##      given factor input (verified directly) - `levels(out) <- levels(mask)`
##      restores it explicitly afterward.
##
##    - `dilate()`: standard label-image dilation - grow existing labels
##      into `NA` cells only, via a majority/mode vote among each `NA`
##      cell's real-valued neighbours (`terra::focal(fun = "modal")`,
##      confirmed to exist and work directly, tie-breaking however
##      terra's own implementation does, not reimplemented here). Cells
##      that already carry a real label are left untouched - `ifel(!is.na(mask),
##      mask, modal_fill)` - this is the genuinely load-bearing design
##      point: unlike binary/continuous dilation (where `max()` naturally
##      never changes an already-`1` cell, since the neighbourhood always
##      includes at least that value), a blind modal filter applied
##      everywhere would ALSO relabel already-classified interior cells
##      near a strong neighbouring majority - a smoothing/denoising
##      operation, not dilation (dilation only ever grows into background,
##      never reassigns existing foreground). Restricting the modal fill
##      to `NA` cells only is what keeps this monotonic, the defining
##      property of dilation.
##
##    Consequence for `tophat()`/`bottomhat()`: the real-valued residual
##    formula above (`mask - opening(mask)`) is meaningless for category
##    codes (subtracting land-cover code 3 from code 1 is not "how much
##    top-hat"). The categorical analogue is a boolean indicator instead -
##    `tophat`: cells that HAD a label but were removed by `opening()`
##    (`!is.na(mask) & is.na(opening(mask))`); `bottomhat`: cells that were
##    `NA` but got FILLED by `closing()` (`is.na(mask) & !is.na(closing(mask))`)
##    - the direct label-image analogue of the binary definition, and
##    identical to it (not just similar) when input actually is binary,
##    since "had a label" / "is 1" and "removed" / "became 0" coincide
##    exactly there. Returned as a plain (non-factor) 0/1 SpatRaster, since
##    the result is inherently a flag, not a category.
##
##    No caveat, no warning: every one of the six functions computes
##    something genuinely well-defined and standard for categorical
##    input, not an arithmetic shortcut applied to arbitrary codes, so
##    there's nothing to warn the caller about.
##
## 3. MISSING structuring-element centre: the anti-/extensivity properties
##    above (`opening(x) <= x`, `closing(x) >= x`, and the categorical
##    erode/dilate reasoning) all assume `kernel` includes its own centre
##    cell - true for every named shortcut at every radius, but a
##    hand-built asymmetric kernel that excludes its own centre can break
##    it. A property of mathematical morphology in general, not specific
##    to this package, and not validated against here.

#' Square structuring element
#'
#' `radius` cells in every direction (so a `2*radius + 1` wide matrix of
#' all `1`s).
#' @param radius integer >= 0
#' @return a plain matrix, all `1`s
#' @export
se_box <- function(radius) {
    matrix(1, 2L * as.integer(radius) + 1L, 2L * as.integer(radius) + 1L)
}

#' Disc (Euclidean ball) structuring element
#'
#' Cell `(dx, dy)` relative to centre is included iff
#' `dx^2 + dy^2 <= radius^2`.
#' @inheritParams se_box
#' @return a `(2*radius+1)`-wide matrix, `1` inside the disc and `NA`
#'   outside it (matching `terra::focal()`'s own convention for a
#'   non-rectangular footprint)
#' @export
se_disc <- function(radius) {
    radius <- as.integer(radius)
    coords <- seq(-radius, radius)
    incl <- outer(coords, coords, function(dx, dy) dx^2 + dy^2 <= radius^2 + 1e-9)
    k <- matrix(as.numeric(incl), nrow = length(coords))
    k[k == 0] <- NA
    k
}

#' Diamond (L1/Manhattan ball) structuring element
#'
#' Cell `(dx, dy)` relative to centre is included iff
#' `abs(dx) + abs(dy) <= radius`.
#' @inheritParams se_box
#' @return a `(2*radius+1)`-wide matrix, `1` inside the diamond and `NA`
#'   outside it
#' @export
se_diamond <- function(radius) {
    radius <- as.integer(radius)
    coords <- seq(-radius, radius)
    incl <- outer(coords, coords, function(dx, dy) abs(dx) + abs(dy) <= radius)
    k <- matrix(as.numeric(incl), nrow = length(coords))
    k[k == 0] <- NA
    k
}

#' @noRd
.mask01 <- function(mask) terra::ifel(is.na(mask), 0, mask)

#' Reject a kernel `terra::focal()` would silently mishandle - see this
#' file's own header for the concrete, verified failure modes (fractional
#' values, values > 1, `0` used instead of `NA`) this closes off. Every
#' non-`NA` entry must be exactly `1`, and at least one entry must be
#' non-`NA` (an all-`NA` kernel has nothing for `focal()` to compute over).
#' @param kernel a matrix
#' @noRd
.check_kernel <- function(kernel) {
    vals <- kernel[!is.na(kernel)]
    if (length(vals) == 0 || sum(kernel, na.rm = TRUE) < 1) {
        stop("`kernel` must have at least one included (`1`) cell - an all-`NA` kernel has nothing ",
             "for the morphological operator to compute over.")
    }
    if (!all(vals == 1)) {
        stop("`kernel` must contain only `1` (include) and `NA` (exclude) values. terra::focal() treats ",
             "any other value as a MULTIPLICATIVE weight, not a boolean include/exclude flag - this ",
             "silently corrupts erode()/dilate()'s output (fractional values, values > 1, or an all-zero ",
             "result if `0` was used in place of `NA`) rather than erroring on its own. Build a custom ",
             "kernel with only `1`/`NA` entries, or use se_box()/se_disc()/se_diamond().")
    }
}

#' Standard label-image erosion: a cell keeps its own category label only
#' if every cell in `kernel`'s footprint shares that exact label,
#' otherwise it becomes `NA`. See this file's own header for the full
#' derivation (equivalent to per-class binary erosion, unioned across
#' classes) and why `na.rm = FALSE` here specifically.
#' @param mask a categorical (`terra::is.factor()`) SpatRaster
#' @param kernel a structuring element matrix (`1`/`NA`)
#' @return categorical SpatRaster, same grid and levels as `mask`
#' @noRd
.erode_core_categorical <- function(mask, kernel) {
    lo <- terra::focal(mask, w = kernel, fun = min, na.rm = FALSE)
    hi <- terra::focal(mask, w = kernel, fun = max, na.rm = FALSE)
    # not `ifel(lo == hi, mask, NA)` - see this file's own header for why
    mismatch <- is.na(lo) | is.na(hi) | (lo != hi)
    out <- terra::ifel(!mismatch, mask, NA)
    # focal()'s min/max strip the factor/levels attribute even on factor
    # input; `levels<-` isn't exported by terra, only registered as an S3
    # method, so this is called unqualified rather than as `terra::levels<-`
    levels(out) <- terra::levels(mask)
    out
}

#' Standard label-image dilation: an already-labelled cell is left
#' untouched; an `NA` cell is filled with the majority label among its
#' real-valued neighbours in `kernel`'s footprint (`NA` if none exist).
#' Restricting the fill to `NA` cells only is what keeps this monotonic -
#' see this file's own header for why a blind modal filter over every
#' cell would be a smoothing operation, not dilation.
#' @inheritParams .erode_core_categorical
#' @return categorical SpatRaster, same grid and levels as `mask`
#' @noRd
.dilate_core_categorical <- function(mask, kernel) {
    filled <- terra::focal(mask, w = kernel, fun = "modal", na.rm = TRUE)
    terra::ifel(!is.na(mask), mask, filled)
}

#' @noRd
.erode_core <- function(mask, kernel) {
    if (terra::is.factor(mask)) return(.erode_core_categorical(mask, kernel))
    # fillvalue = 0, not focal()'s own NA default - see this file's own
    # header (rspatial/terra#243) for why; dilate() needs no equivalent fix
    terra::focal(.mask01(mask), w = kernel, fun = min, na.rm = TRUE, fillvalue = 0)
}

#' @noRd
.dilate_core <- function(mask, kernel) {
    if (terra::is.factor(mask)) return(.dilate_core_categorical(mask, kernel))
    terra::focal(.mask01(mask), w = kernel, fun = max, na.rm = TRUE)
}

#' Morphological erosion of a terra SpatRaster
#'
#' For binary (0/1) or continuous input, the local minimum over `kernel`'s
#' footprint at every cell (grayscale erosion, a strict generalization of
#' binary erosion, not a different operation); for a categorical
#' (`terra::is.factor()`) input, standard label-image erosion - a cell
#' keeps its own label only if its entire neighbourhood shares it,
#' otherwise it becomes `NA` (see this file's own header for the full
#' reasoning behind both).
#' @param mask binary, continuous, or categorical SpatRaster
#' @param kernel a structuring element: any matrix, or one of `se_box()`,
#'   `se_disc()`, `se_diamond()`. `1` includes a cell in the footprint,
#'   `NA` excludes it (matching `terra::focal()`'s own `w` convention) -
#'   any other value errors (see `.check_kernel()`'s own comments for why).
#' @return SpatRaster, same grid as `mask` and the same value type
#'   (binary/continuous/categorical) as `mask` itself
#' @examples
#' r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, crs = "local")
#' terra::values(r) <- 0
#' r[8:13, 8:13] <- 1
#' erode(r, kernel = se_disc(1))
#' @export
erode <- function(mask, kernel = se_box(1)) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_kernel(kernel)
    .erode_core(mask, kernel)
}

#' Morphological dilation of a terra SpatRaster
#'
#' For binary (0/1) or continuous input, the local maximum over `kernel`'s
#' footprint at every cell (grayscale dilation, a strict generalization of
#' binary dilation, not a different operation); for a categorical
#' (`terra::is.factor()`) input, standard label-image dilation -
#' already-labelled cells are left untouched, and `NA` cells are filled
#' with the majority label among their real-valued neighbours (see this
#' file's own header for the full reasoning behind both).
#' @inheritParams erode
#' @return SpatRaster, same grid and value type as `mask` (see `erode()`)
#' @examples
#' r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, crs = "local")
#' terra::values(r) <- 0
#' r[8:13, 8:13] <- 1
#' dilate(r, kernel = se_disc(1))
#' @export
dilate <- function(mask, kernel = se_box(1)) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_kernel(kernel)
    .dilate_core(mask, kernel)
}

#' Morphological opening (erode then dilate) of a terra SpatRaster
#'
#' Removes small bright features and thin protrusions without shrinking
#' the overall shape the way `erode()` alone does. Generalizes to
#' continuous input as grayscale opening, and to categorical input as
#' label-image opening (removes small/thin single-class regions), the
#' same operation in each case, via `erode()`/`dilate()`'s own dispatch.
#' @inheritParams erode
#' @return SpatRaster, same grid and value type as `mask` (see `erode()`)
#' @export
opening <- function(mask, kernel = se_box(1)) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_kernel(kernel)
    .dilate_core(.erode_core(mask, kernel), kernel)
}

#' Morphological closing (dilate then erode) of a terra SpatRaster
#'
#' Fills small dark gaps/notches without growing the overall shape the way
#' `dilate()` alone does. Generalizes to continuous input as grayscale
#' closing, and to categorical input as label-image closing (fills small
#' gaps within a class region), the same operation in each case, via
#' `erode()`/`dilate()`'s own dispatch.
#' @inheritParams erode
#' @return SpatRaster, same grid and value type as `mask` (see `erode()`)
#' @export
closing <- function(mask, kernel = se_box(1)) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_kernel(kernel)
    .erode_core(.dilate_core(mask, kernel), kernel)
}

#' Top-hat (white-hat) transform of a terra SpatRaster
#'
#' Highlights small bright features/protrusions narrower than the
#' structuring element. For binary/continuous input, the real-valued
#' residual `mask - opening(mask)` (reduces exactly to the classic
#' `mask AND NOT opening(mask)` when input is 0/1). For categorical input,
#' a boolean indicator instead - `1` at cells that had a label but lost it
#' under `opening()` - since subtracting category codes is meaningless
#' (see this file's own header for why).
#'
#' Binary/continuous input specifically: `NA` outside `mask`'s own
#' footprint, not a numeric `0` - `opening()`'s own output never actually
#' contains `NA` for non-categorical input (`.mask01()` replaces `NA`
#' with `0` internally so `focal()` treats it as background), so without
#' this restriction the residual would be a fully-defined `0` arbitrarily
#' far outside the shape. Restoring `NA` there loses no information:
#' `mask - opening(mask)` is exactly `0` everywhere outside `mask`'s own
#' footprint regardless (`opening(x) <= x` pointwise, and opening never
#' creates value beyond `mask`'s own extent).
#' @inheritParams erode
#' @return SpatRaster, same grid as `mask`. Binary/continuous input: `NA`
#'   outside `mask`'s own footprint, non-negative and real-valued inside
#'   it wherever `kernel` includes its own centre (true for
#'   `se_box()`/`se_disc()`/`se_diamond()` at every radius - a custom
#'   kernel that excludes its own centre can violate this, a property of
#'   grayscale morphology in general). Categorical input: plain `0`/`1`,
#'   not itself categorical - the result is a flag, not a category.
#' @export
tophat <- function(mask, kernel = se_box(1)) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_kernel(kernel)
    o <- .dilate_core(.erode_core(mask, kernel), kernel)
    if (terra::is.factor(mask)) return(terra::ifel(!is.na(mask) & is.na(o), 1, 0))
    terra::ifel(is.na(mask), NA, .mask01(mask) - o)
}

#' Bottom-hat (black-hat) transform of a terra SpatRaster
#'
#' Highlights small dark gaps/notches narrower than the structuring
#' element. For binary/continuous input, the real-valued residual
#' `closing(mask) - mask` (reduces exactly to the classic
#' `closing(mask) AND NOT mask` when input is 0/1). For categorical input,
#' a boolean indicator instead - `1` at cells that were unlabelled (`NA`)
#' but got filled by `closing()` - since subtracting category codes is
#' meaningless (see this file's own header for why).
#'
#' Binary/continuous input specifically: `NA` only where `mask` was
#' ITSELF `NA` (genuinely missing data, not just background) AND beyond
#' one `kernel` radius of `mask`'s own footprint - narrower than
#' `tophat()`'s equivalent restriction, deliberately, on both counts.
#' Not "wherever `mask` is `NA`" the way `tophat()` is: `bottomhat()`'s
#' entire purpose is flagging small gaps ADJACENT to (not inside) the
#' shape, and a gap IS an `NA` cell in `mask` - restoring `NA` there would
#' erase exactly the cells this function exists to highlight. And not
#' "wherever `mask` is confirmed `0`" either: a plain binary (`0`/`1`, no
#' `NA` at all) mask has no genuinely missing data anywhere, so its
#' `bottomhat()` output stays `0`/`1` everywhere, matching the classic
#' definition exactly. The one-kernel-radius reach bound is exact, not a
#' heuristic margin: `closing()` can only ever differ from background
#' within one kernel radius of `mask`'s own footprint (dilate grows by at
#' most one radius, the subsequent erode can only give at most that much
#' back).
#' @inheritParams erode
#' @return SpatRaster, same grid as `mask` (see `tophat()`'s own note on
#'   value type and non-negativity; `NA` region described above, not the
#'   same one)
#' @export
bottomhat <- function(mask, kernel = se_box(1)) {
    on.exit(.cleanup_tmpfiles(), add = TRUE)
    .check_kernel(kernel)
    c_ <- .erode_core(.dilate_core(mask, kernel), kernel)
    if (terra::is.factor(mask)) return(terra::ifel(is.na(mask) & !is.na(c_), 1, 0))
    valid <- !is.na(mask) & (mask != 0)
    reach <- .dilate_core(terra::ifel(valid, 1, 0), kernel) > 0
    terra::ifel(is.na(mask) & !reach, NA, c_ - .mask01(mask))
}
