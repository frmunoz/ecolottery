\name{coalesc_abc}
\alias{coalesc_abc}
\alias{do.simul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimation of neutral and non-neutral parameters of community assembly using
  Approximate Bayesian Computation (ABC)
}
\description{
  Estimates parameters of neutral migration-drift dynamics (through migration
  rate m and parameters of environmental filtering (through a filtering function
  \code{filt.abc()}) from the composition of a local community and the related
  regional pool. 
}
\usage{
  coalesc_abc(comm.obs, pool = NULL, multi = "single", traits = NULL,
              f.sumstats, filt.abc = NULL, params = NULL, 
              theta.max = NULL, nb.samp = 10^6, parallel = TRUE,
              tol = NULL, pkg = NULL, method="rejection")
  do.simul(J, pool = NULL, multi = "single", nb.com = NULL, 
              traits = NULL, f.sumstats = NULL, filt.abc = NULL,
              params, theta.max = NULL, nb.samp = 10^6, 
              parallel = TRUE, tol = NULL, pkg = NULL,
              method = "rejection")
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
  \item{pkg}{
    packages needed for calculation of \code{filt.abc} and/or \code{f.sumstats}.
  }
  \item{method}{
    the method to be used in ABC estimation (see help on
    \code{\link[=abc]{abc()}} function  of package \code{abc} for further
    information). Can be "rejection" (default), "loclinear", "neuralnet" or "ridge".
  }
  \item{J}{
    local community size.
  }
  \item{nb.com}{
    number of communities.
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

}
\author{
  F. Munoz
}
\examples{
# Trait-dependent filtering function
filt_gaussian <- function(t, params) exp(-(t-params[1])^2/(2*params[2]^2))

# Definition of parameters and their range
params <- data.frame(rbind(c(0, 1), c(0.05, 1)))
row.names(params) <- c("topt", "sigmaopt")
# Number of values to sample in prior distributions
nb.samp <- 10^6 # Should be large

\dontrun{
# Basic summary statistics
f.sumstats <- function(com) array(dimnames=list(c("cwm", "cwv", "cws",
                                                  "cwk", "S", "Es")),
                                  c(mean(com[,3]), var(com[,3]), 
                                    e1071::skewness(com[,3]),  
                                    e1071::kurtosis(com[,3]),
                                    vegan::specnumber(table(com[,2])),
                                    vegan::diversity(table(com[,2]))))

# An observed community is here simulated (known parameters)
comm <- coalesc(J = 400, m = 0.5, theta = 50,
                filt = function(x) filt_gaussian(x, c(0.2, 0.1)))

# ABC estimation of the parameters based on observed community composition 
## Warning: this function may take a while
res <- coalesc_abc(comm$com, comm$pool, f.sumstats = f.sumstats,
                   filt.abc = filt_gaussian, params = params, 
                   nb.samp = nb.samp, parallel = TRUE, 
                   pkg = c("e1071","vegan"), method = "neuralnet")
plot(res$abc, param = res$par)
hist(res$abc)

# Cross validation
## Warning: this function is slow
res$cv <- abc::cv4abc(param = res$par, sumstat = res$ss, nval = 1000,
                      tols = c(0.01, 0.1, 1), method = "neuralnet")
plot(res$cv)

# Multiple community option
# When the input is a site-species matrix, use argument multi="tab"
# See vignette Barro_Colorado for more details

# When the input is a list of communities, use argument multi="seqcom"
comm.obs <- list()

comm.obs[[1]] <- cbind(rep(1,400), coalesc(J = 400, m = 0.5, filt = function(x)
                                           filt_gaussian(x, c(0.2, 0.1)),
                                           pool = comm$pool)$com))
comm.obs[[2]] <- cbind(rep(2,400), coalesc(J = 400, m = 0.5, filt = function(x) 
                                           filt_gaussian(x, c(0.5, 0.1)),
                                           pool = comm$pool)$com))
comm.obs[[3]] <- cbind(rep(3,400), coalesc(J = 400, m = 0.5, filt = function(x) 
                                           filt_gaussian(x, c(0.8, 0.1)),
                                           pool = comm$pool)$com))

comm.obs <- lapply(comm.obs, as.matrix)

res <- coalesc_abc(comm.obs, comm$pool, multi="seqcom", f.sumstats=f.sumstats,
                   filt.abc = filt_gaussian, params = params, nb.samp = nb.samp,
                   parallel = TRUE, pkg = c("e1071","vegan"), tol = 0.1,
                   method = "neuralnet")
                   
lapply(res$abc, summary)

}
}

\keyword{coalescent}
\keyword{Approximate Bayesian Computation}
\keyword{niche-based dynamics}
\keyword{neutral dynamics}
