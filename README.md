# ecolottery
[![Build Status](https://travis-ci.org/frmunoz/ecolottery.svg?branch=master)](https://travis-ci.org/frmunoz/ecolottery)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/GITHUB_USERNAME/REPO?branch=master&svg=true)](https://ci.appveyor.com/project/frmunoz/ecolottery)

Coalescent-Based Simulation of Ecological Communities

This package allows quick and efficient simulation of communities with combination of environmental filtering and neutral birth-death dynamics. We also provide an Approximate Bayesian Computation tool to assess parameters of these processes from observed communities.
We propose a coalescent-based approach that presents significant advantages over alternative algorithms, allowing intensive calculation schemes. It does not require simulating community dynamics from an initial state forward in time, and the genealogical approach can incorporate the influence of environmental filtering on immigration and establishment. The sampling process allows simulating varying sample size in a flexible way.
The package also includes a forward-in-time simulation algorithm to investigate a broader range of niche-based dynamics, albeit it requires greater resources and computation time.

## Installing the package

The version 1.0.0 is available on CRAN,
https://cran.r-project.org/web/packages/ecolottery/index.html
```r
install.packages("ecolottery")
```

The development version on Github can be installed using the `devtools` package:
```r
devtools::install_github("frmunoz/ecolottery/pkg", build_vignettes = TRUE)
```

## How to use the package

The two main functions of the package are `coalesc()` and `forward()`. You can learn more on how to use them by looking at the [introduction vignette](pkg/vignettes/coalesc_vignette.Rmd).

If you want to see how `coalesc_abc()` allows estimating parameters of community assembly with the ABC approach, have a look at the [Barro Colorado vignette](pkg/vignettes/Barro_Colorado.Rmd).
