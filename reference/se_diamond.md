# Diamond (L1/Manhattan ball) structuring element

Cell `(dx, dy)` relative to centre is included iff
`abs(dx) + abs(dy) <= radius`.

## Usage

``` r
se_diamond(radius)
```

## Arguments

- radius:

  integer \>= 0

## Value

a `(2*radius+1)`-wide matrix, `1` inside the diamond and `NA` outside it
