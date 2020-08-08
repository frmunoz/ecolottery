\name{coalesc_easyABC}
\alias{coalesc_easyABC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimation of parameters of community assembly using
  Approximate Bayesian Computation (ABC), with algorithms of optimal exploration of prior distributions
}
\description{
  Estimates parameters of neutral migration-drift dynamics (through migration
  rate m and parameters of environmental filtering (through a filtering function
  \code{filt.abc()}) from the composition of a local community and the related
  regional pool. 
}
\usage{
 coalesc_easyABC(comm.obs, pool = NULL, multi = "single", prop = F, traits = NULL,
  f.sumstats, filt.abc = NULL, filt.vect = F, migr.abc = NULL, size.abc = NULL, add = F, 
  var.add = NULL, params = NULL, par.filt = NULL, par.migr = NULL, 
  par.size = NULL, constr = NULL, scale = F, dim.pca = NULL, svd = F, 
  theta.max = NULL, nb.samp = 10^6, parallel = TRUE, nb.core = NULL, 
  tol = NULL, pkg = NULL, type = NULL, method.seq = "Lenormand", 
  method.mcmc = "Marjoram_original", method.abc = "rejection", alpha = 0.5)
}
\arguments{
  \item{comm.obs}{
    the observed community composition. If \code{multi = "single"} (default),
    should be a matrix or data.frame of individuals on rows with their
    individual id (first column), and species id (second column).
  }
  \item{pool}{
    composition of the regional pool to which the local community is
    hypothesized to be related through migration dynamics with possible
    environmental filtering. Should be a matrix of individuals on rows with
    their individual id (first column), species id (second column), and
    (optionally) the trait values of the individuals.
    When \code{multi = "tab"} or \code{"seqcom"}, different pools can be 
    set for each community: \code{pool} must be then defined as a list of 
    matrices with the same columns as above.
  }
  \item{multi}{
    structure of the community inputs:
    \itemize{
      \item if \code{multi = "single"}, \code{comm.obs} represents a single
      community
      
      \item if \code{multi = "tab"}, the user provides a site-species matrix
      (\strong{sites in rows and species in columns})
      
      \item if \code{multi = "seqcom"}, \code{comm.obs} contains a list of
      communities. If \code{prop = F}, the local community composition is given
      with one individual per row (first column is individual id, second column is 
      species id, subsequent columns provide trait values).
      If \code{prop = F}, the local community composition is given with one species
      per row (first column if species id, second column is coverage, 
      subsequent columns provide trait values).
    }
  }
  \item{prop}{
    indicates if the community composition is given in term of relative species abundances, 
    cover or proportions. In this case, a parameter of effective community size is estimated.
    Default is FALSE. If no prior is provided for \code{J} based on argumeng \code{par.size}, 
    it is by default a uniform  distribution between 100 and 1000.
  }
  \item{traits}{
    the trait values of species in the regional pool. It is used if trait
    information is not provided in \code{pool}. In this case, intraspecific
    trait variation is assumed to be null. Species names of \code{pool} must
    be included in \code{traits}.
  }
  \item{f.sumstats}{
    a function calculating the summary statistics of community
    composition. It is be used to compare observed and simulated community
    composition in ABC estimation. The first input argument of the function concerns community
    composition, under a format depending on \code{multi} and on \code{prop}. 
    If \code{multi="comm"} this first argument is the local community composition with
    individual and species labels as first columns and (optionally) trait values as subsequent
    columns. 
    There can be a second argument taking the average trait values per species (derived from
    \code{pool} if not provided in \code{traits}, see above), and \code{var.add} 
    (see below) can be sent as a third argument. See examples of different options below.
    The function should return a vector of summary statistics.
  }
  \item{filt.abc}{
    the hypothesized environmental filtering function. It is a function of
    individual trait values (first argument), additional parameters to be estimated and (optionally) 
    environmental variables defined in \code{var.add}.
    If \code{NULL}, neutral communities will be simulated and only m will be estimated. 
  }
   \item{filt.vect}{
    indicates whether the filtering function can be vectorized. It means that the function
    can take as input a vector of trait values and provide a vector of the corresponding 
    weights.
  }
  \item{migr.abc}{
    a function defining immigration probability in local communities. It can be a function of 
    environmental variables defined in \code{var.add} and of parameters to be estimated.
  }
  \item{size.abc}{
    a function defining local community sizes. It can be used when \code{prop = T}.
    It can be a function of environmental variables defined in \code{var.add} and of 
    parameters to be estimated.
  }
   \item{add}{
    indicates if additional variables must be passed to \code{filt.abc} and/or \code{migr.abc}
    . It can be, for instance, environmental data conditioning the trait-based filtering 
    and/or immigration in the community. Default is FALSE.
  }
    \item{var.add}{
    additional variables to be passed to \code{filt.abc} and/or \code{migr.abc} 
    when \code{add = T}.
  }
  \item{params}{
    equivalent to par.filt (see below): it is kept in the function for compatibility with 
    previous versions.
  }
  \item{par.filt}{
    a matrix of the bounds of the parameters used in \code{filt.abc}. The row
    names of \code{par.filt} provide the parameter names used in ABC calculation
    and output. First column contains minimum values and second column contains
    maximum values.
  }
  \item{par.migr}{
    a matrix of the bounds of the parameters used in \code{migr.abc}. The row
    names of \code{par.migr} provide the parameter names used in ABC calculation
    and output. First column contains minimum values and second column contains
    maximum values.
  }
  \item{par.size}{
    a matrix of the bounds of the parameters used in \code{size.abc}. The row
    names of \code{par.size} provide the parameter names used in ABC calculation
    and output. First column contains minimum values and second column contains
    maximum values.
  }
  \item{constr}{
    constraints that must be set on parameter values, i.e., relationships that must
    be met when drawing the values in prior distributions. It must be a vector of 
    character strings each including an operator relating some of the parameters.
    The names of parameters must be consistent with those used in \code{par.filt},
    \code{par.migr} and \code{par.size}. Default is \code{NULL} (no constraint).
  }
  \item{scale}{
    should the summary statistics be scaled. Default is \code{FALSE}.
  }
  \item{dim.pca}{
    gives a number of independent dimensions to calculate from the simulated summary statistics, 
    by performing Principal Component Analysis (PCA). Default is \code{NULL}, in which case 
    no PCA is performed. 
    If \code{svd = T}, Singular Value Decomposition (SVD) will be performed instead of PCA.
  }
  \item{svd}{
    indicates if Singular Value Decomposition (SVD) must be performed to get the basic 
    dimensions requested in \code{dim.pca}. Will be ignored if \code{dim.pca = NULL}.
  }
  \item{theta.max}{
    if \code{pool = NULL}, regional abundances will be simulated following a
    log-series distribution. The function will estimate the \code{theta}
    parameter of this distribution. \code{theta.max} then provides the upper
    bound for this estimation.
  }
  \item{nb.samp}{
   the number of parameter values to be sampled in ABC calculation. Random
   values of parameters of environmental filtering (see \code{filt.abc} and
   \code{params}) and of migration (denoted as \code{m}) are drawn from a
   uniform distribution between minimum and maximum values provided in
   \code{params} (and between 0 and 1 for m).
  }
  \item{parallel}{
   boolean. If \code{parallel = TRUE}, the function will perform parallel
   processing using the \code{\link[=parLapply]{parLapply()}} function of
   package \code{parallel}. Default is \code{parallel = FALSE}.
  }
  \item{nb.core}{
   number of cores to be used in parallel computation if \code{parallel = TRUE}. 
   If \code{NULL} (default), all the cores minus 1 are used.
   If \code{1}, only one core is used (i.e., no parallel computing).
  }
  \item{tol}{
    the tolerance value used in ABC estimation (see help in
    \code{\link[=abc]{abc()}} function of package \code{abc} for further
    information). Default is \code{NULL}.
  }
  \item{pkg}{
    packages needed for calculation of \code{filt.abc} and/or \code{f.sumstats}.
  }
  \item{type}{
    the type of algorithm to be used in EasyABC. Can be either "standard" (using package abc, default)" "seq" (sequential), "mcmc" or "annealing". Three later options are based on the EasyABC package.
  }
  \item{method.seq}{
    when type = "seq", gives the algorithm for sequential sampling scheme, which is passed to \code{\link[=EasyABC]{
    ABC_sequential}}. Can be "Lenormand" (Default), "Drovandi", "Delmoral", "Beaumont" or "Emulation".
  }
  \item{method.mcmc}{
    when type = "mcmc", gives the algorithm for MCMC sampling scheme, which is passed to \code{\link[=EasyABC]{
    ABC_mcmc}}. Can be "Marjoram_original"(Default), "Marjoram" or "Wegmann". The method "Marjoram_original" cannot      be used with multiple cores.
  }
 \item{method.abc}{
    the method to be used in ABC estimation (see help on
    \code{\link[=abc]{abc()}} function  of package \code{abc} for further
    information). Can be "rejection", "loclinear", "neuralnet" or "ridge". Default is "rejection".
  }
  \item{alpha}{
    a positive number between 0 and 1 (strictly) used when performing sequential ABC method. \code{alpha} is the proportion of particles rejected at each step in the algorithm "Drovandi". This is the proportion of particles kept at each step in the algorithms "Delmoral", "Lenormand" and "Emulation". Default value is 0.5
  }

}
\details{
\code{coalesc_easyABC()} provides ABC estimation of assembly parameters for one (if \code{multi = "single"},
default) or several communities (if \code{multi = "tab" or "seqcom"}) related to the same
regional pool. The function uses optimization algorithms implemented in the package \code{EasyABC}.
}
\value{
  \item{par}{
    parameter values used in simulations.
  }
  \item{obs}{
    observed summary statistics.
  }
  \item{obs.scaled}{
    observed summary statistics standardized according to the mean and standard
    deviation of simulated values.
  }
  \item{ss}{
    standardized summary statistics of the communities simulated with parameter
    values listed in \code{par}.
  }
  \item{ss.scale}{
    data frame including the mean and the standard deviation used for standardization of  
    observed and summary statistics.
  }
  \item{abc}{
    a single (if \code{multi = FALSE}, default) or a list of abc objects
    including ABC estimation information for each community provided in input
    (\code{comm.obs}). It is of class \code{abc} for \code{type="standard"} and of class \code{EasyABC} otherwise.
  }
}
\seealso{
 \code{\link[=coalesc_abc]{coalesc_abc()}}
}
\references{
  Jabot, F., T. Faure, and N. Dumoulin 2013. EasyABC: performing efficient
  approximate Bayesian computation sampling schemes using R. Methods in 
  Ecology and Evolution 4:684-687.
}
\author{
  F. Munoz, E. Barthelemy
}

\keyword{coalescent}
\keyword{Approximate Bayesian Computation}
\keyword{environmental filtering}
\keyword{neutral dynamics}
\keyword{EasyABC}
