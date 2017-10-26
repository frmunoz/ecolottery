\name{coalesc_abc2}
\alias{coalesc_abc2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimation of neutral and non-neutral parameters of community assembly using
  Approximate Bayesian Computation (ABC), using the algorithms from the EasyABC package
}
\description{
  Estimates parameters of neutral migration-drift dynamics (through migration
  rate m and parameters of environmental filtering (through a filtering function
  \code{filt.abc()}) from the composition of a local community and the related
  regional pool. 
  Compared to \code{coalesc_abc}, \code{coalesc_abc2} uses the algorithms implemented
  in the package \code{EasyABC} for ABC processing.
}
\usage{
  coalesc_abc2(comm.obs, pool, multi = "single", traits = NULL, f.sumstats, filt.abc, 
  params,  theta.max = NULL, nb.samp = 10^6, parallel = T, tol = 10^-4, type = "seq", 
  method.seq = "Beaumont", method.mcmc = "Marjoram",
  method.abc = "neuralnet")
}
\arguments{
  \item{comm.obs}{
    the observed community composition. If \code{multi = FALSE} (default),
    should be a matrix or data.frame of individuals on rows with their
    individual id (first column), and species id (second column).
  }
  \item{pool}{
    composition of the regional pool to which the local community is
    hypothesized to be related through migration dynamics with possible
    environmental filtering. Should be a matrix of individuals on rows with
    their individual id (first column), species id (second column), and
    (optionally) the trait values of the individuals.
  }
  \item{multi}{
    structure of the community inputs:
    \itemize{
      \item if \code{multi = "single"}, \code{comm.obs} contains a single
      community
      
      \item if \code{multi = "tab"}, the user provides a site-species matrix
      (\strong{sites in rows and species in columns})
      
      \item if \code{multi = "seqcom"}, \code{comm.obs} contains a list of
      communities
    }
  }
  \item{traits}{
    the trait values of species in the regional pool. It is used if trait
    information is not provided in \code{pool}. In this case, intraspecific
    trait variation is assumed to be null.
  }
  \item{f.sumstats}{
    a function allowing to calculate the summary statistics of local community
    composition. Will be used to compare observed and simulated community
    composition in the ABC estimation. It should take a community as input and
    output a list of summary statistics.
  }
  \item{filt.abc}{
    the hypothesized environmental filtering function. It is a function of
    individual trait values and additional parameters to be estimated.
  }
  \item{params}{
    a matrix of the bounds of the parameters used in \code{filt.abc}. The row
    names of \code{params} provide the parameter names used in ABC calculation
    and output. First column contains minimum values and second column contains
    maximum values.
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
   package \code{parallel}.
  }
  \item{tol}{
    the tolerance value used in ABC estimation (see help in
    \code{\link[=abc]{abc()}} function of package \code{abc} for further
    information).
  }
  %\item{pkg}{
  %  packages needed for calculation of \code{filt.abc} and/or \code{f.sumstats}.
  %}
  %\item{J}{
  %  local community size.
  %}
  %\item{nb.com}{
  %  number of communities.
  %}
  \item{type}{
    the type of algorithm to be used in EasyABC. Can be either "seq" (sequential, default), "mcmc" or "annealing".
  }
  \item{method.seq}{
    when type = "seq", gives the algorithm for sequential sampling scheme, which is passed to \code{\link[=EasyABC]{
    ABC_sequential}}. Can be "Beaumont" (Default), "Drovandi", "Delmoral", "Lenormand" or "Emulation".
  }
  \item{method.mcmc}{
    when type = "mcmc", gives the algorithm for MCMC sampling scheme, which is passed to \code{\link[=EasyABC]{
    ABC_mcmc}}. Can be "Marjoram_original", "Marjoram" (Default) or "Wegmann". The method "Marjoram_original" cannot      be used with multiple cores.
  }
 \item{method.abc}{
    the method to be used in ABC estimation (see help on
    \code{\link[=abc]{abc()}} function  of package \code{abc} for further
    information). Default is "neuralnet".
  }

}
\details{
\code{coalesc_abc()} performs ABC estimation for one (if \code{multi = FALSE},
default) or several communities (if \code{multi = TRUE}) related to the same
regional pool.

\code{do.simul()} provides the simulated communities used in ABC estimation,
and is not intended to be used directly.

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
    data frame including the mean and the standard deviation used for standardization of       observed and summary statistics.
  }
  \item{abc}{
    a single (if \code{multi = FALSE}, default) or a list of abc objects
    including ABC estimation information for each community provided in input
    (\code{comm.obs}).
  }
}
\seealso{
 \code{\link[=abc]{abc()}} in \code{abc} package,
 \code{\link[=parLapply]{parLapply()}} in \code{parallel} package.
}
\references{
  Jabot, F., and J. Chave. 2009. Inferring the parameters of the neutral theory
  of biodiversity using phylogenetic information and implications for tropical
  forests. Ecology Letters 12:239-248.

  Csillery, K., M. G. B. Blum, O. E. Gaggiotti, and O. Francois. 2010.
  Approximate Bayesian computation (ABC) in practice. Trends in Ecology &
  Evolution 25:410-418.

  Csillery, K., O. Francois, and M. G. Blum. 2012. abc: an R package for
  Approximate Bayesian Computation (ABC). Methods in Ecology and Evolution
  3:475-479.
  
  Jabot, F., T. Faure, and N. Dumoulin 2013. EasyABC: performing efficient
  approximate Bayesian computation sampling schemes using R. Methods in 
  Ecology and Evolution 4:684-687.
}
\author{
  F. Munoz
}
\examples{
# Simulated community
J <- 500
m <- 0.1
# Define a regional pool of species with equal abundances
pool <- cbind(1:10000, rep(1:500, 20), rep(NA, 10000))
# Uniform distribution of trait values
t.sp <- runif(500)
# No intraspecific variation
pool[,3] <- t.sp[pool[,2]]
# Function for environmental filtering
sigma <- 0.1
filt_gaussian <- function(x, par) exp(-(x - par[1])^2/(2*par[2]^2))
# Stabilizing habitat filtering around t = 0.1
comm <- coalesc(J, m, filt = function(x) filt_gaussian(x, c(0.1,sigma)), pool = pool)

\dontrun{
f.sumstats <- function(com) array(dimnames=list(c("cwm", "cwv", "cws",
                                                  "cwk", "S", "Es")),
                                  c(mean(com[,3]), var(com[,3]), 
                                    e1071::skewness(com[,3]),  
                                    e1071::kurtosis(com[,3]),
                                    vegan::specnumber(table(com[,2])),
                                    vegan::diversity(table(com[,2]))))
f.sumstats(comm$com)
# Definition of parameters and their range
params <- data.frame(rbind(c(0, 1), c(0.05, 1)))
row.names(params) <- c("topt", "sigmaopt")
# Number of values to sample in prior distributions
nb.samp <- 1000 # Should be large

comm.abc2 <- coalesc_abc2(comm$com, comm$pool, f.sumstats=f.sumstats, filt.abc=filt_gaussian, params=params, nb.samp=nb.samp)
}
}

\keyword{coalescent}
\keyword{Approximate Bayesian Computation}
\keyword{niche-based dynamics}
\keyword{neutral dynamics}
\keyword{EasyABC}
