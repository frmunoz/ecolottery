# lottery
[![Build Status](https://travis-ci.org/frmunoz/lottery.svg?branch=master)](https://travis-ci.org/frmunoz/lottery)

Coalescent-Based Simulation of Ecological Communities

This package allows quick and efficient simulation of communities undergoing both neutral and niche-based assembly processes. We also provide an Approximate Bayesian Computation tool to assess parameters of these processes from observed communities.
We propose a coalescent-based approach that presents significant advantages over alternative algorithms, allowing intensive calculation schemes. It does not require simulating community dynamics from an initial state forward in time, and the genealogical approach can incorporate the influence of environmental filtering and species niche differences. The sampling process allows simulating varying sample size in a flexible way.


## Installing the package

As the package is not yet on CRAN you can install the Github version using the `devtools` package:
```r
devtools::install_github("frmunoz/lottery/pkg", build_vignettes = TRUE)
```

## How to use the package

The two main functions of the package are `coalesc()` and `forward()`. You can learn more on how to use them by looking at the [introduction vignette](pkg/vignettes/coalesc_vignette.Rmd).

If you want to see how `coalesc_abc()` allows estimating parameters of community assembly with the ABC approach, have a look at the [Barro Colorado vignette](pkg/vignettes/Barro_Colorado.Rmd).
