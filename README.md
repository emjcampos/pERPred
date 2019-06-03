
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pERPred

<!-- badges: start -->

<!-- badges: end -->

The goal of pERPred is to implement the algorithm proposed in “Principle
ERP reduction and analysis: Estimating and using principle ERP waveforms
underlying ERPs across subjects, channels, and conditions” authored by
Emilie Campos and Chad Hazlett et al (link to paper). The algorithm
finds a set of basis functions upon which ERP analysis can be applied.

## Installation

You can install the released version of pERPred from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("pERPred")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("emjcampos/pERPred")
```

## Example

An example analysis would begin with a dataframe structured with one
column for Task, Subject, and Time, and a column for each electrode.

``` r
library(pERPred)
```

``` r
simulated_data
#> # A tibble: 37,500 x 28
#>     Task Subject  Time electrode_1 electrode_2 electrode_3 electrode_4
#>    <int>   <int> <dbl>       <dbl>       <dbl>       <dbl>       <dbl>
#>  1     1       1 0.004   -1.12           0.846       0.139       0.106
#>  2     1       1 0.008    0.912         -0.524       2.14       -1.13 
#>  3     1       1 0.012   -1.05           0.897       1.93       -1.01 
#>  4     1       1 0.016   -1.81           0.460       0.348      -1.37 
#>  5     1       1 0.02    -0.650          0.324      -1.58       -1.37 
#>  6     1       1 0.024   -0.905          1.44        1.76       -1.21 
#>  7     1       1 0.028   -2.59           1.26        1.08       -1.06 
#>  8     1       1 0.032   -1.32          -0.915      -0.822      -1.40 
#>  9     1       1 0.036    0.486         -0.930       0.455      -1.28 
#> 10     1       1 0.04    -0.000448      -0.251      -1.79       -1.36 
#> # … with 37,490 more rows, and 21 more variables: electrode_5 <dbl>,
#> #   electrode_6 <dbl>, electrode_7 <dbl>, electrode_8 <dbl>,
#> #   electrode_9 <dbl>, electrode_10 <dbl>, electrode_11 <dbl>,
#> #   electrode_12 <dbl>, electrode_13 <dbl>, electrode_14 <dbl>,
#> #   electrode_15 <dbl>, electrode_16 <dbl>, electrode_17 <dbl>,
#> #   electrode_18 <dbl>, electrode_19 <dbl>, electrode_20 <dbl>,
#> #   electrode_21 <dbl>, electrode_22 <dbl>, electrode_23 <dbl>,
#> #   electrode_24 <dbl>, electrode_25 <dbl>
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
