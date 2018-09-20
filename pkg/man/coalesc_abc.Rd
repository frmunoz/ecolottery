\name{coalesc_abc}
\alias{coalesc_abc}
\alias{do.simul.coalesc}
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
  coalesc_abc(comm.obs, pool = NULL, multi = "single", prop = F, traits = NULL,
              f.sumstats, filt.abc = NULL, migr.abc = NULL, size.abc = NULL, add = F, 
              var.add = NULL, params = NULL, par.filt = NULL, par.migr = NULL, 
              par.size = NULL, constr = NULL, scale = T, dim.pca = NULL, svd = F, 
              theta.max = NULL, nb.samp = 10^6, parallel = TRUE, nb.core = NULL, 
              tol = NULL, pkg = NULL, method = "rejection")
  do.simul.coalesc(J, pool = NULL, multi = "single", prop = F, nb.com = NULL, 
              traits = NULL, f.sumstats = NULL, filt.abc = NULL, migr.abc = NULL, 
              size.abc = NULL, add = F, var.add = NULL, params = NULL, par.filt = NULL, 
              par.migr = NULL, par.size = NULL, constr = NULL, dim.pca = NULL, svd = F, 
              theta.max = NULL, nb.samp = 10^6, parallel = TRUE, nb.core = NULL, 
              tol = NULL, pkg = NULL, method = "rejection")
}
\arguments{
  \item{comm.obs}{
    the observed community composition. If \code{multi = "single"} (default),
    should be a matrix or data.frame of individuals on rows with their
    individual id (first column), and species id (second column).
  }
  \item{pool}{
    give the composition of the regional pool to which a local community is
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
      \item if \code{multi = "single"}, \code{comm.obs} contains a single
      community
      
      \item if \code{multi = "tab"}, the user provides a site-species matrix
      (\strong{sites in rows and species in columns})
      
      \item if \code{multi = "seqcom"}, \code{comm.obs} contains a list of
      communities
    }
  }
  \item{prop}{
    indicates if the community composition is given in term of relative species abundances, 
    cover or proportions. In this case, a parameter of effective community size is estimated.
    Presently only available with multi = "tab". Default is FALSE. If no prior is provided for
    \code{J}, it is by default a uniform distribution between 100 and 1000.
  }
  \item{traits}{
    the trait values of species in the regional pool. It is used if trait
    information is not provided in \code{pool}. In this case, intraspecific
    trait variation is assumed to be null. Species names of \code{pool} must
    be included in \code{traits}.
  }
  \item{f.sumstats}{
    a function allowing to calculate the summary statistics of local community
    composition. Will be used to compare observed and simulated community
    composition in the ABC estimation. It should take a community as input and
    output a list of summary statistics.
  }
  \item{filt.abc}{
    the hypothesized environmental filtering function. It is a function of
    individual trait values, additional parameters to be estimated and (optionally) 
    environmental variables defined in \code{var.add}.
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
    be met when drawing the values in prior distributions. It must be a list of 
    character strings each including an operator relating some of the parameters.
    The names of parameters must be consistent with those used in \code{par.filt},
    \code{par.migr} and \code{par.size}. Default is \code{NULL} (no constraint).
  }
  \item{scale}{
    should the summary statistics be scaled. Default is \code{TRUE}.
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
   package \code{parallel}.
  }
  \item{nb.core}{
   number of cores to be used in parallel computation if \code{parallel = TRUE}. 
   If \code{NULL}, all the cores minus 1 are used.
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

\code{do.simul.coalesc()} provides the simulated communities used in ABC estimation,
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
    data frame including the mean and the standard deviation used for standardization of
    observed and summary statistics.
  }
  \item{abc}{
    a single (if \code{multi = FALSE}, default) or a list of abc objects
    including ABC estimation information for each community provided in input
    (\code{comm.obs}).
  }
  \item{call}{
    the call function.
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
\dontrun{

# 1/ Analysis of community assembly in a single community

# Trait-dependent filtering function
filt_gaussian <- function(t, par) exp(-(t - par[1])^2/(2*par[2]^2))

# Definition of parameters of this function and their range
par.range <- data.frame(rbind(c(0, 1), c(0.05, 1)))
row.names(par.range) <- c("topt", "sigmaopt")

# Basic summary statistics describing local community composition
f.sumstats <- function(com) array(dimnames=list(c("cwm", "cwv", "cws",
                                                  "cwk", "S", "Es")),
                                  c(mean(com[,3]), var(com[,3]), 
                                    e1071::skewness(com[,3]),  
                                    e1071::kurtosis(com[,3]),
                                    vegan::specnumber(table(com[,2])),
                                    vegan::diversity(table(com[,2]))))

# An observed community is simulated with known parameter values
comm <- coalesc(J = 400, m = 0.5, theta = 50,
                filt = function(x) filt_gaussian(x, c(0.2, 0.1)))

# ABC estimation of the parameters of filtering and migration based on observed 
# community composition 
# Number of values to sample in prior distributions
nb.samp <- 1000 # Should be large
## Warning: this function may take a while
res.single <- coalesc_abc(comm$com, comm$pool, f.sumstats = f.sumstats,
                   filt.abc = filt_gaussian, par.filt = par.range, 
                   nb.samp = nb.samp, parallel = TRUE, tol = 0.5,
                   pkg = c("e1071","vegan"), method = "neuralnet")
plot(res.single$abc, param = res.single$par)
hist(res.single$abc)

# Cross validation
## Warning: this function is slow
res.single$cv <- abc::cv4abc(param = res.single$par, sumstat = res.single$ss, nval = 100,
                      tols = c(0.01, 0.1, 1), method = "neuralnet")
plot(res.single$cv)

# 2/ Analysis of community assembly in multiple communities with a common pool of 
# immigrants

# When the input is a site-species matrix, use argument multi="tab"

# 2.1/ Environment-dependent environmental filtering

# The variation of optimal trait values can depend on environmental variation across communities
filt_gaussian_env <- function(t, par, env) exp(-(t-par[1]*env+par[2])^2/(2*par[3]^2))
# Vector of environmental conditions across communities
env1 <- matrix(runif(20))

# Simulate a set of 20 communities with varying environmental filtering in
# different environmental conditions
# A common source pool of immigrants
pool <- coalesc(J = 10000, theta = 50)$pool
meta1 <- c()
for(i in 1:20)
  meta1[[i]] <- coalesc(J = 400, pool = pool, m = 0.5, 
                filt = function(x) filt_gaussian_env(x, c(0.5, 0.2, 0.1), env1[i,]))
# Build a species-by-site table
tab1 <- array(0, c(20, max(pool$sp)))
for(i in 1:20) {
  compo <- table(meta1[[i]]$com$sp)
  tab1[i, as.numeric(names(compo))] <- compo
}
colnames(tab1) <- as.character(1:ncol(tab1))
tab1 <- tab1[, colSums(tab1)!=0]

# Definition of parameters and their range
# We will estimate the slope a and intercept b of the linear relationship between 
# optimal trait values and environmental conditions
par.range <- data.frame(rbind(c(-1, 1), c(0, 1), c(0.05, 1)))
row.names(par.range) <- c("a", "b", "sigmaopt")

# Basic summary statistics
# The function will provide trait-based statistics and taxonomic diversity in
# each community
f.sumstats <- function(tab, traits) array(dimnames=
  list(c(paste(rep("cwm",nrow(tab)),1:nrow(tab)), paste(rep("cwv",nrow(tab)),1:nrow(tab)), 
  paste(rep("cws",nrow(tab)),1:nrow(tab)), paste(rep("cwk",nrow(tab)),1:nrow(tab)), 
  paste(rep("S",nrow(tab)),1:nrow(tab)), paste(rep("Es",nrow(tab)),1:nrow(tab)), 
  paste(rep("beta",nrow(tab)),1:nrow(tab)))),
  c(apply(tab, 1, function(x) Weighted.Desc.Stat::w.mean(traits[colnames(tab),], x)), 
  apply(tab, 1, function(x) Weighted.Desc.Stat::w.var(traits[colnames(tab),], x)), 
  apply(tab, 1, function(x) Weighted.Desc.Stat::w.skewness(traits[colnames(tab),], x)),  
  apply(tab, 1, function(x) Weighted.Desc.Stat::w.kurtosis(traits[colnames(tab),], x)),
  apply(tab, 1, function(x) sum(x!=0)),
  apply(tab, 1, vegan::diversity),
  colMeans(as.matrix(vegan::betadiver(tab,"w")))))

# ABC estimation of the parameters of trait-environment relationship and of migration,
# based on observed community composition 
## Warning: this function may take a while
res.tab1 <- coalesc_abc(tab1, pool, multi = "tab", f.sumstats = f.sumstats,
                   filt.abc = filt_gaussian_env, par.filt = par.range, 
                   add = T, var.add = env1,
                   nb.samp = nb.samp, parallel = TRUE, tol = 0.5,
                   pkg = c("e1071","vegan"), method = "neuralnet")

plot(res.tab1$abc, param=res.tab1$par)                   
hist(res.tab1$abc)

# 2.2/ Environment-dependent environmental filtering and immigration

# In this case, the migration rate depends of another environmental variable
# representing, e.g., community isolation
# There will be two environmental variables used for parameter inference, one conditioning
# environmental filtering, and the other conditioning migration
env2 <-  matrix(runif(20))
env <- cbind(env1, env2)
colnames(env) <- c("env1", "env2")
# Migration depends on environmental conditions following some linear relationship
migr_env <- function(par, env) (par[2]-par[1])*env+par[1]

# Simulate a set of 20 communities with environment-dependent environmental filtering 
# and immigration
meta2 <- c()
for(i in 1:20)
  meta2[[i]] <- coalesc(J = 400, pool = pool, m = migr_env(c(0.25,0.75), env[i,2]), 
                filt = function(x) filt_gaussian_env(x, c(0.5, 0.2, 0.1), env[i,1]))
# Build a species-by-site table
tab2 <- array(0, c(20, max(pool$sp)))
for(i in 1:20) {
  compo <- table(meta2[[i]]$com$sp)
  tab2[i, as.numeric(names(compo))] <- compo
}
colnames(tab2) <- as.character(1:ncol(tab2))
tab2 <- tab2[, colSums(tab2)!=0]

# Definition of parameters and their range
# We will estimate the slope a and intercept b of the linear relationship between 
# optimal trait values and environmental variable 1, and slope c and intercept d in the
# logistic function determining the variation of migration rate with environmental variable 2
par.filt.range <- data.frame(rbind(c(-1, 1), c(0, 1), c(0.05, 1)))
row.names(par.filt.range) <- c("a", "b", "sigmaopt")
par.migr.range <- data.frame(rbind(c(0, 1), c(0, 1)))
row.names(par.migr.range) <- c("c", "d")

# ABC estimation of the parameters of trait-environment and migration-environment 
# relationships, based on observed community composition 
## Warning: this function may take a while
nb.samp <- 200
res.tab2 <- coalesc_abc(tab2, pool, multi = "tab", f.sumstats = f.sumstats,
                   filt.abc = function(x, par, env) filt_gaussian_env(x, par, env[1]), 
                   migr.abc = function(par, env) migr_env(par, env[2]),
                   par.filt = par.filt.range, par.migr = par.migr.range,
                   add = T, var.add = env,
                   nb.samp = nb.samp, parallel = FALSE, tol = 0.5,
                   pkg = c("e1071","vegan"), method = "neuralnet")

plot(res.tab2$abc, param=res.tab2$par)                   
hist(res.tab2$abc)

# Check result consistency
comm <- 1
coeff <- summary(res.tab2$abc)
plot(cbind(sapply(pool$tra, function(x) filt_gaussian_env(x, c(0.5, 0.2, 0.1), env[comm,1])), 
sapply(pool$tra, function(x) filt_gaussian_env(x, coeff[3,1:3], env[comm,1]))), 
xlab = paste("Expected filter value in community ", comm, sep=" "), ylab = "Estimated filter value")
abline(0,1)
plot(cbind(sapply(env[,2], function(x) migr_env(c(0.25, 0.75), x)), 
sapply(env[,2], function(x) migr_env(coeff[3,4:5], x))), 
xlab = "Expected migration value", ylab = "Estimated migration value")
abline(0,1)

# 3/ Distinct pools of migrants across communities

# The immigrants can be drawn from distinct pools for each community, to represent a local context.
# In this case, the pool argument sent to coalesc_abc must be a list of the local pools.

# Simulate several communities with distinct pools and same environmental filter
pool.loc <- c()
meta3 <- c()
theta <- 2.5*(1:20) # The pools differ in diversity
# We first define a baseline pool of immigrants common to all communities
baseline <- coalesc(J = 10000, theta = 25)$pool
baseline$sp <- paste(0, baseline$sp, sep = "_")
for(i in 1:20)
{
  pool.loc[[i]] <- coalesc(J = 10000, theta = theta[i])$pool
  pool.loc[[i]]$sp <- paste(i, pool.loc[[i]]$sp, sep = "_")
  pool.loc[[i]] <- rbind(pool.loc[[i]], baseline)
  meta3[[i]] <- coalesc(J = 400, pool = pool.loc[[i]], m = 0.5, 
                filt = function(x) filt_gaussian(x, c(0.25, 0.1)))
}
# Build a species-by-site table
tab3 <- array(0, c(20, length(unique(Reduce(rbind, pool.loc)[,2]))))
colnames(tab3) <- unique(Reduce(rbind, pool.loc)[,2])
for(i in 1:20) {
  compo <- table(meta3[[i]]$com$sp)
  tab3[i, names(compo)] <- compo
}
tab3 <- tab3[, colSums(tab3)!=0]

# ABC estimation of the parameters of trait-environment and migration-environment 
# relationships, based on observed community composition 
## Warning: this function may take a while
par.range <- data.frame(rbind(c(0, 1), c(0.05, 1)))
row.names(par.range) <- c("topt", "sigmaopt")
res.tab3 <- coalesc_abc(tab3, pool.loc, multi = "tab", f.sumstats = f.sumstats,
                   filt.abc = filt_gaussian,
                   par.filt = par.range, 
                   nb.samp = nb.samp, parallel = FALSE, tol = 0.5,
                   pkg = c("e1071","vegan"), method = "neuralnet")

plot(res.tab3$abc, param=res.tab3$par)  
hist(res.tab3$abc)
}
}

\keyword{coalescent}
\keyword{Approximate Bayesian Computation}
\keyword{niche-based dynamics}
\keyword{neutral dynamics}
