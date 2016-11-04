\name{coalesc_abc}
\alias{coalesc_abc}
\alias{do.simul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimation of neutral and non-neutral parameters of community assembly using Approximate Bayesian Computation (ABC)
}
\description{
  Estimates parameters of neutral migration-drift dynamics (through migration rate m and parameters of environmental filtering (through a filtering function filt.abc) from the composition of a local community and the related regional pool. 
}
\usage{
  coalesc_abc(comm.obs, pool, multi = F, traits = NULL, f.sumstats, filt.abc, params, nb.samp = 10^6, parallel = T, tol = 1*10^-4,
              pkg = NULL, method="neuralnet")
  do.simul(J, pool, traits = NULL, f.sumstats, filt.abc, params, nb.samp = 10^6, parallel = T, tol = 1*10^-4, pkg = NULL, 
           method = "neuralnet")
}
\arguments{
  \item{comm.obs}{
    The observed community composition. If multi = F (default), should be a matrix of individuals on rows with their individual id (first column), and species id (second column).
  }
  \item{pool}{
    Composition of the regional pool to which the local community is hypothesized to be related through migration dynamics with possible environmental filtering. Should be a matrix of individuals on rows with their individual id (first column), species id (second column), and (optionnally) the trait values of the individuals.
  }
  \item{multi}{
  If TRUE, the user provides in comm.obs multiple community data, i.e., the first column of comm.obs includes community number, the second column the individual id and the third the species id.
  }
  \item{traits}{
    The trait values of species in the regional pool. It is used if trait information is not provided in pool. In this case, intraspecific trait variation is assumed to be null.
  }
  \item{f.sumstats}{
    A function allowing to calculate the summary statistics of local community composition. Will be used to compare observed and simulated community composition in the ABC estimation.
  }
  \item{filt.abc}{
    The hypothesized environmental filtering function. It is a function of individual trait values and additional parameters to be estimated.
  }
  \item{params}{
    A matrix of the bounds of the parameters used in filt.abc. The row names of params provide the parameter names used in ABC calculation and output.
  }
  \item{nb.samp}{
    The number of parameter values to be sampled in ABC calculation. Random values of parameters of environmental filtering (see filt.abc and params) and of migration (denoted as m) are drawn from a uniform distribution between minimum and maximum values provided in params (and between 0 and 1 for m).
  }
  \item{parallel}{
    If TRUE, the function will perform parallel processing using the parLapply function of package parallel.
  }
  \item{tol}{
    The tolerance value used in ABC estimation (see help on abc function for further information).
  }
  \item{pkg}{
    The packages that need to be loaded for calculation of filt.abc and/or f.sumstats.
  }
  \item{method}{
    The method to be used in ABC estimation (see help on abc function for further information).
  }
  \item{J}{
    Local community size.
  }
}
\details{
coalesc_abc performes ABC estimation for one (if multi = F, default) or several communities (if multi = T) related to the same regional pool.

do.simul provides the simulated communities used in ABC estimation, and is not intended to be used directly.

}
\value{
  \item{par}{
    The parameter values used in simulations.
  }
  \item{par}{
    Summary statistics of the communities simulated with parameter values listed in par.
  }
  \item{abc}{
    A single (if multi = F) or a list of abc objects including ABC estimation information for each community provided in input (comm.obs).
  }
}
\references{
  Jabot, F., and J. Chave. 2009. Inferring the parameters of the neutral theory of biodiversity using phylogenetic information and implications for tropical forests. Ecology Letters 12:239-248.

  Csilléry, K., M. G. B. Blum, O. E. Gaggiotti, and O. François. 2010. Approximate Bayesian computation (ABC) in practice. Trends in Ecology & Evolution 25:410-418.

  Csilléry, K., O. François, and M. G. Blum. 2012. abc: an R package for Approximate Bayesian Computation (ABC). Methods in Ecology and Evolution 3:475–479.

}
\author{
  F. Munoz
}

\examples{
# Trait-dependent filtering function
filt_gaussian <- function(t,params) exp(-(t-params[1])^2/(2*params[2]^2))

# Definition of parameters and their range
params <- data.frame(rbind(c(0,1),c(0.05,1)))
row.names(params) <- c("topt","sigmaopt")
# Number of values to sample in prior distributions
nb.samp <- 10^6; # Should be large

## Not run:
# Basic summary statistics
f.sumstats <- function(com,t) array(dimnames=list(c("cwm","cwv","cws","cwk","S","Es")),c(mean(t),
      var(t),
      e1071::skewness(t),  
      e1071::kurtosis(t),
      vegan::specnumber(table(com)),
      vegan::diversity(table(com))))

# An observed community is here simulated (known parameters)
comm <- coalesc(J = 400, m = 0.5, theta = 50, filt = function(x) filt_gaussian(x,c(0.2,0.1)))

# ABC estimation of the parameters based on observed community composition
res <- coalesc_abc(comm$com, comm$pool, f.sumstats=f.sumstats, filt.abc=filt_gaussian, 
                   params=params, nb.samp = nb.samp, parallel = T, pkg =c("e1071","vegan"))
plot(res$abc,param=res$par)
hist(res$abc)

# Multiple community option
comm.obs <- c();
comm.obs <- rbind(comm.obs,cbind(rep(1,400),as.matrix(coalesc(J = 400, m = 0.5, filt = function(x) 
  filt_gaussian(x,c(0.2,0.1)),pool=comm$pool)$com)))
comm.obs <- rbind(comm.obs,cbind(rep(2,400),as.matrix(coalesc(J = 400, m = 0.5, filt = function(x) 
  filt_gaussian(x,c(0.5,0.1)),pool=comm$pool)$com)))
comm.obs <- rbind(comm.obs,cbind(rep(3,400),as.matrix(coalesc(J = 400, m = 0.5, filt = function(x) 
  filt_gaussian(x,c(0.8,0.1)),pool=comm$pool)$com)))
res <- coalesc_abc(comm.obs, comm$pool, multi=T, f.sumstats=f.sumstats, filt.abc=filt_gaussian,
                   params=params, nb.samp = nb.samp, parallel = T, pkg =c("e1071","vegan"))
summary(res$abc[[1]])
summary(res$abc[[2]])
summary(res$abc[[3]])
## End(Not run)
}

\keyword{coalescent}
\keyword{Approximate Bayesian Computation}
\keyword{niche-based dynamics}
\keyword{neutral dynamics}
