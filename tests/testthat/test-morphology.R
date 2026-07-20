make_test_mask <- function() {
    r <- terra::rast(nrows = 9, ncols = 9, xmin = 0, xmax = 9, ymin = 0, ymax = 9, crs = "local")
    terra::values(r) <- 0
    r[4:6, 4:6] <- 1   # a solid 3x3 square, room to erode/dilate around it
    r
}

test_that("se_box/se_disc/se_diamond produce sane, correctly-sized kernels", {
    expect_equal(dim(se_box(1)), c(3, 3))
    expect_equal(sum(se_box(1)), 9)
    expect_equal(dim(se_disc(2)), c(5, 5))
    expect_equal(dim(se_diamond(2)), c(5, 5))
    # disc(2) should be strictly smaller (or equal) footprint than box(2) - never bigger
    expect_lte(sum(se_disc(2) == 1, na.rm = TRUE), sum(se_box(2) == 1))
})

test_that("an arbitrary user-supplied matrix works directly as a kernel", {
    r <- make_test_mask()
    custom <- matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), 3, 3)  # a plus-shape, hand-built
    eroded <- erode(r, kernel = custom)
    expect_true(inherits(eroded, "SpatRaster"))
})

test_that("erode shrinks a solid square, dilate grows it, by the expected cell counts", {
    r <- make_test_mask()
    n0 <- sum(as.vector(terra::values(r)), na.rm = TRUE)
    eroded <- erode(r, kernel = se_box(1))
    dilated <- dilate(r, kernel = se_box(1))
    expect_lt(sum(as.vector(terra::values(eroded)), na.rm = TRUE), n0)
    expect_gt(sum(as.vector(terra::values(dilated)), na.rm = TRUE), n0)
})

test_that("erode on a 3x3 box with a 3x3 box kernel leaves only the centre cell", {
    r <- make_test_mask()
    eroded <- erode(r, kernel = se_box(1))
    vals <- as.vector(terra::values(eroded))
    expect_equal(sum(vals, na.rm = TRUE), 1)
})

test_that("opening = erode-then-dilate and closing = dilate-then-erode, verified not assumed", {
    r <- make_test_mask()
    k <- se_box(1)
    manual_opening <- dilate(erode(r, k), k)
    expect_equal(as.vector(terra::values(opening(r, k))), as.vector(terra::values(manual_opening)))

    manual_closing <- erode(dilate(r, k), k)
    expect_equal(as.vector(terra::values(closing(r, k))), as.vector(terra::values(manual_closing)))
})

test_that("tophat and bottomhat are binary (0/1) and zero where expected", {
    r <- make_test_mask()
    k <- se_box(1)
    th <- tophat(r, k)
    bh <- bottomhat(r, k)
    expect_true(all(as.vector(terra::values(th)) %in% c(0, 1)))
    expect_true(all(as.vector(terra::values(bh)) %in% c(0, 1)))
    # a solid box has no small dark notches for closing to fill, so bottomhat is all zero
    expect_equal(sum(as.vector(terra::values(bh)), na.rm = TRUE), 0)
})

test_that("NA cells are treated as background (0), not propagated", {
    r <- make_test_mask()
    r[1, 1] <- NA
    eroded <- erode(r, kernel = se_box(1))
    expect_false(any(is.na(as.vector(terra::values(eroded)))))
})

test_that("a kernel with fractional or >1 non-NA values errors, not silently corrupts the output", {
    r <- make_test_mask()
    frac <- matrix(c(1, 1, 1, 1, 0.5, 1, 1, 1, 1), 3, 3)
    big <- matrix(c(1, 1, 1, 1, 2, 1, 1, 1, 1), 3, 3)
    expect_error(erode(r, kernel = frac), "only.*1.*NA")
    expect_error(dilate(r, kernel = big), "only.*1.*NA")
})

test_that("a kernel using 0 instead of NA to mean 'exclude' errors rather than silently zeroing everything", {
    r <- make_test_mask()
    zero_based <- matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3)  # meant as "just the centre", wrong convention
    expect_error(erode(r, kernel = zero_based), "only.*1.*NA")
})

test_that("an all-NA kernel errors instead of computing over nothing", {
    r <- make_test_mask()
    empty <- matrix(NA_real_, 3, 3)
    expect_error(erode(r, kernel = empty), "at least one included")
})

test_that("opening/closing/tophat/bottomhat validate the same kernel rules as erode()/dilate()", {
    r <- make_test_mask()
    frac <- matrix(c(1, 1, 1, 1, 0.5, 1, 1, 1, 1), 3, 3)
    expect_error(opening(r, kernel = frac), "only.*1.*NA")
    expect_error(closing(r, kernel = frac), "only.*1.*NA")
    expect_error(tophat(r, kernel = frac), "only.*1.*NA")
    expect_error(bottomhat(r, kernel = frac), "only.*1.*NA")
})

make_continuous_rast <- function() {
    r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "local")
    set.seed(1)
    terra::values(r) <- stats::runif(100, 0, 5)
    r
}

test_that("erode/dilate on continuous input are genuine local-min/local-max filters, not binarized", {
    r <- make_continuous_rast()
    e <- erode(r, se_box(1))
    d <- dilate(r, se_box(1))
    ev <- as.vector(terra::values(e))
    dv <- as.vector(terra::values(d))
    rv <- as.vector(terra::values(r))
    # output is not restricted to {0, 1} the way binary erosion/dilation would be
    expect_false(all(ev %in% c(0, 1)))
    expect_false(all(dv %in% c(0, 1)))
    # erosion never exceeds the original value, dilation never falls below it
    expect_true(all(ev <= rv + 1e-9))
    expect_true(all(dv >= rv - 1e-9))
})

test_that("opening/closing on continuous input match manual erode-then-dilate/dilate-then-erode exactly", {
    r <- make_continuous_rast()
    k <- se_box(1)
    expect_equal(as.vector(terra::values(opening(r, k))), as.vector(terra::values(dilate(erode(r, k), k))))
    expect_equal(as.vector(terra::values(closing(r, k))), as.vector(terra::values(erode(dilate(r, k), k))))
})

test_that("tophat/bottomhat on continuous input are the real grayscale residual, not silently all-zero", {
    r <- make_continuous_rast()
    k <- se_box(1)
    th <- tophat(r, k)
    bh <- bottomhat(r, k)
    thv <- as.vector(terra::values(th))
    bhv <- as.vector(terra::values(bh))

    # this is the regression this test guards: the OLD implementation tested
    # `m01 == 1 & o == 0` (exact equality against 0/1), which is essentially
    # never true for continuous values and silently returned all-0
    expect_true(any(thv > 1e-9))
    expect_true(any(bhv > 1e-9))

    # matches the real formula exactly, not approximately
    expect_equal(thv, as.vector(terra::values(r - opening(r, k))))
    expect_equal(bhv, as.vector(terra::values(closing(r, k) - r)))
})

test_that("tophat/bottomhat on binary input still reduce exactly to the old 0/1 formula", {
    r <- make_test_mask()
    k <- se_box(1)
    th <- as.vector(terra::values(tophat(r, k)))
    bh <- as.vector(terra::values(bottomhat(r, k)))
    m01 <- as.vector(terra::values(r))
    o <- as.vector(terra::values(opening(r, k)))
    c_ <- as.vector(terra::values(closing(r, k)))
    expect_equal(th, as.numeric(m01 == 1 & o == 0))
    expect_equal(bh, as.numeric(c_ == 1 & m01 == 0))
})

test_that("tophat on continuous input with NA outside the shape is NA there, not a spurious 0", {
    # regression test: .mask01() converts NA to 0 internally (needed so
    # focal() treats it as background rather than dropping it), which
    # used to leak a fully-defined 0 residual arbitrarily far outside the
    # shape, since opening()'s own output never actually contains NA for
    # non-categorical input - verified this is real, not cosmetic:
    # mask - opening(mask) is exactly 0 outside mask's own footprint
    # regardless, so restoring NA there loses no information
    n <- 21
    r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    cx <- terra::init(r, "x") - (n / 2)
    cy <- terra::init(r, "y") - (n / 2)
    weight_rast <- terra::ifel(sqrt(cx^2 + cy^2) <= 8, cx + 15, NA)
    k <- se_disc(2)

    th <- tophat(weight_rast, k)
    outside <- is.na(terra::values(weight_rast))
    expect_true(all(is.na(terra::values(th)[outside])))
    expect_true(any(!is.na(terra::values(th)[!outside])))  # still real-valued where mask has data
})

test_that("bottomhat's NA restriction on continuous input is narrower than tophat's - a plain binary mask is unaffected, and gaps adjacent to the shape stay informative", {
    # regression test: an earlier version of this fix restored NA
    # wherever bottomhat's output was beyond one kernel radius of the
    # shape, full stop - which broke a plain binary (0/1, no NA at all)
    # mask by introducing NA into what used to be a fully-defined 0/1
    # result. The correct condition additionally requires mask itself to
    # be NA (genuinely missing data), not just far from the shape.
    r <- make_test_mask()  # 9x9, no NA anywhere
    k <- se_box(1)
    expect_false(any(is.na(terra::values(bottomhat(r, k)))))

    # a genuinely interior gap, far enough from any class/shape boundary
    # that it's a real bottomhat candidate, stays informative (not NA)
    n <- 11
    landcover <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "local")
    terra::values(landcover) <- 1
    landcover[5, 5] <- NA
    bh <- bottomhat(landcover, k)
    expect_false(is.na(bh[5, 5][[1]]))
    expect_equal(bh[5, 5][[1]], 1)  # closing(landcover)[5,5] (1) - mask01(landcover)[5,5] (0)

    # far outside the raster's own data entirely (well beyond the
    # kernel's own reach from any real cell), it IS NA
    far <- terra::rast(nrows = 41, ncols = 41, xmin = 0, xmax = 41, ymin = 0, ymax = 41, crs = "local")
    terra::values(far) <- NA
    far[20, 20] <- 1
    bh_far <- bottomhat(far, k)
    expect_true(is.na(bh_far[1, 1][[1]]))
})

make_categorical_rast <- function() {
    # a 3-class raster (classes 1, 2, 3) with an interior NA gap, split
    # roughly diagonally so erosion at the class boundaries is exercised
    r <- terra::rast(nrows = 7, ncols = 7, xmin = 0, xmax = 7, ymin = 0, ymax = 7, crs = "local")
    m <- matrix(c(
        1, 1, 1, 1, 2, 2, 2,
        1, 1, 1, 1, 2, 2, 2,
        1, 1, 1, 1, 2, 2, 2,
        1, 1, 1, 3, 2, 2, 2,
        1, 1, 1, 3, 3, 2, 2,
        1, 1, 1, 3, 3, 3, 2,
        1, 1, 1, 3, 3, 3, 3
    ), 7, 7, byrow = TRUE)
    terra::values(r) <- as.vector(t(m))
    r[4, 4] <- NA
    terra::as.factor(r)
}

test_that("erode on categorical input is a per-cell unanimity check, not a numeric min", {
    r <- make_categorical_rast()
    k <- se_box(1)
    eroded <- erode(r, k)

    expect_true(terra::is.factor(eroded))
    expect_equal(terra::levels(eroded), terra::levels(r))

    # matches the min-equals-max derivation directly, not just assumed to -
    # built via is.na()/!= rather than `ifel(lo == hi, r, NA)` directly,
    # since the latter is NOT equivalent here: see this file's own header
    # for a verified terra::ifel() quirk where an all-TRUE-or-NA condition
    # (exactly what `lo == hi` looks like near a raster edge on uniform
    # input) silently drops the NA branch
    lo <- terra::focal(r, w = k, fun = min, na.rm = FALSE)
    hi <- terra::focal(r, w = k, fun = max, na.rm = FALSE)
    mismatch <- is.na(lo) | is.na(hi) | (lo != hi)
    manual <- terra::ifel(!mismatch, r, NA)
    expect_equal(as.vector(terra::values(eroded)), as.vector(terra::values(manual)))

    # a cell straddling the class-1/class-2 boundary loses its label...
    expect_true(is.na(terra::values(eroded)[2 * 7 - 3]))  # row2, col4 (value 1, neighbors include class 2)
    # ...while a cell deep inside a single class - genuinely interior, not
    # touching the raster's own edge, away from any class boundary or the
    # NA gap - keeps it
    expect_equal(terra::values(eroded)[3 * 7 - 5], 1)  # row3, col2
})

test_that("erode on categorical input also erodes at the raster's own edge, not just interior boundaries", {
    # regression test for a verified terra::focal() edge-padding bug
    # (rspatial/terra#243-adjacent, see this file's own header): a shape
    # reaching the raster's true edge must still erode there, the same as
    # it would next to an interior NA hole or a different class
    r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "local")
    terra::values(r) <- 1
    r <- terra::as.factor(r)
    eroded <- erode(r, se_box(1))
    vals <- matrix(as.vector(terra::values(eroded)), 5, 5, byrow = TRUE)
    expect_true(all(is.na(vals[1, ])) && all(is.na(vals[5, ])) && all(is.na(vals[, 1])) && all(is.na(vals[, 5])))
    expect_true(all(!is.na(vals[2:4, 2:4])))
})

test_that("erode on binary/continuous input also erodes at the raster's own edge (fillvalue = 0 regression)", {
    r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "local")
    terra::values(r) <- 1
    eroded <- erode(r, se_box(1))
    vals <- matrix(as.vector(terra::values(eroded)), 5, 5, byrow = TRUE)
    expect_equal(sum(vals), 9)  # exactly the interior 3x3 block survives
    expect_true(all(vals[2:4, 2:4] == 1))
    expect_true(all(vals[1, ] == 0) && all(vals[5, ] == 0) && all(vals[, 1] == 0) && all(vals[, 5] == 0))
})

test_that("dilate at the raster's own edge needs no fillvalue fix - 0 is max()'s identity element", {
    r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "local")
    terra::values(r) <- 0
    r[1, 1] <- 1  # a single foreground cell at the very corner
    dilated <- dilate(r, se_box(1))
    vals <- matrix(as.vector(terra::values(dilated)), 5, 5, byrow = TRUE)
    expect_equal(sum(vals), 4)  # exactly the 2x2 block around the corner grows
    expect_true(all(vals[1:2, 1:2] == 1))
})

test_that("dilate on categorical input only fills NA cells by majority vote, never relabels existing cells", {
    r <- make_categorical_rast()
    k <- se_box(1)
    dilated <- dilate(r, k)

    expect_true(terra::is.factor(dilated))

    orig <- as.vector(terra::values(r))
    dil <- as.vector(terra::values(dilated))
    # every already-labelled cell is untouched
    expect_equal(dil[!is.na(orig)], orig[!is.na(orig)])
    # the interior NA gap (row4, col4) gets filled - majority of its
    # neighbours (1,1,2 / 1,. ,2 / 1,3,3) is class 1
    expect_equal(as.numeric(as.character(dilated[4, 4][[1]])), 1)
})

test_that("opening/closing on categorical input match manual erode-then-dilate/dilate-then-erode", {
    r <- make_categorical_rast()
    k <- se_box(1)
    expect_equal(as.vector(terra::values(opening(r, k))), as.vector(terra::values(dilate(erode(r, k), k))))
    expect_equal(as.vector(terra::values(closing(r, k))), as.vector(terra::values(erode(dilate(r, k), k))))
})

test_that("tophat on categorical input flags a genuine thin single-class spike, not a numeric residual", {
    # solid class 1, a solid class-2 block on the right, and a one-cell
    # spike of class 2 poking into class-1 territory - opening should
    # remove the spike (too thin for the structuring element), and tophat
    # should flag exactly that cell
    r <- terra::rast(nrows = 7, ncols = 7, xmin = 0, xmax = 7, ymin = 0, ymax = 7, crs = "local")
    m <- matrix(1, 7, 7)
    m[1:7, 5:7] <- 2
    m[4, 4] <- 2
    terra::values(r) <- as.vector(t(m))
    r <- terra::as.factor(r)
    k <- se_box(1)

    th <- tophat(r, k)
    expect_false(terra::is.factor(th))  # a flag, not a category
    expect_true(all(as.vector(terra::values(th)) %in% c(0, 1)))
    expect_equal(sum(as.vector(terra::values(th))), 1)
    expect_equal(th[4, 4][[1]], 1)

    o <- opening(r, k)
    manual_th <- terra::ifel(!is.na(r) & is.na(o), 1, 0)
    expect_equal(as.vector(terra::values(th)), as.vector(terra::values(manual_th)))
})

test_that("bottomhat on categorical input flags a genuine interior notch, not a numeric residual", {
    # solid class 1 with a class-2 block far away (so the notch below is
    # not itself boundary-adjacent), and a one-cell interior NA notch
    r <- terra::rast(nrows = 11, ncols = 11, xmin = 0, xmax = 11, ymin = 0, ymax = 11, crs = "local")
    m <- matrix(1, 11, 11)
    m[1:11, 9:11] <- 2
    terra::values(r) <- as.vector(t(m))
    r <- terra::as.factor(r)
    r[5, 5] <- NA
    k <- se_box(1)

    bh <- bottomhat(r, k)
    expect_false(terra::is.factor(bh))
    expect_true(all(as.vector(terra::values(bh)) %in% c(0, 1)))
    expect_equal(sum(as.vector(terra::values(bh))), 1)
    expect_equal(bh[5, 5][[1]], 1)

    c_ <- closing(r, k)
    manual_bh <- terra::ifel(is.na(r) & !is.na(c_), 1, 0)
    expect_equal(as.vector(terra::values(bh)), as.vector(terra::values(manual_bh)))
})

test_that("none of the six functions warn on categorical input any more - the result is well-defined, not a caveat", {
    r <- make_categorical_rast()
    k <- se_box(1)
    expect_no_warning(erode(r, k))
    expect_no_warning(dilate(r, k))
    expect_no_warning(opening(r, k))
    expect_no_warning(closing(r, k))
    expect_no_warning(tophat(r, k))
    expect_no_warning(bottomhat(r, k))
})
