# Disc (Euclidean ball) structuring element

Cell `(dx, dy)` relative to centre is included iff
`dx^2 + dy^2 <= radius^2`.

## Usage

``` r
se_disc(radius)
```

## Arguments

- radius:

  integer \>= 0

## Value

a `(2*radius+1)`-wide matrix, `1` inside the disc and `NA` outside it
(matching
[`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)'s
own convention for a non-rectangular footprint)
