# ecolottery - Backward- and forward-in-time simulation of ecological communities
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ecolottery)](https://cran.r-project.org/package=ecolottery)
![](http://cranlogs.r-pkg.org/badges/grand-total/ecolottery?color=yellowgreen)
[![Build Status](https://travis-ci.org/frmunoz/ecolottery.svg?branch=master)](https://travis-ci.org/frmunoz/ecolottery)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/frmunoz/ecolottery?branch=master&svg=true)](https://ci.appveyor.com/project/frmunoz/ecolottery)

The package allows quick and efficient simulation of community composition under diverse combinations of environmental filtering and neutral birth-death dynamics. The package offers a coalescent-based, backward-in-time algorithm, named *coalesc*.
The package also includes a forward-in-time simulation algorithm, named *forward*, to investigate a broader set of niche-based dynamics, albeit at the expense of greater resources and computation time.

*coalesc_abc* is another core function of the package, to perform Approximate Bayesian Computation (ABC) and infer parameters of the processes from observed community composition. This function relies on the *coalesc* function for intensive simulation of biodiversity parrtterns in ABC analyses.

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
