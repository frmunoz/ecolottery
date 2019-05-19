\name{coalesc}
\alias{coalesc}
% -
\title{
Coalescent-based simulation of ecological communities undergoing both neutral and niche-base dynamics
}
\description{
Simulates the composition of a community based on immigration from a regional pool, habitat filtering depending on local environment and species traits, and local birth-death stochastic dynamics. 
}
\usage{
coalesc(J, m = 1, theta = NULL, filt = NULL, add = FALSE,  var.add =NULL, pool = NULL, 
        traits = NULL, Jpool = 50 * J, verbose = FALSE)
}
\arguments{
 \item{J}{
number of individuals in the local community.
}
 \item{m}{
migration rate (if \code{m = 1} the community is a subsample of the regional pool).
}
 \item{theta}{
parameter of neutral dynamics in the regional pool (used only if \code{pool=NULL}), it is the \dQuote{fundamental biodiversity number} (\eqn{\theta}).
}
 \item{filt}{
a function representing the effect of local habitat filtering. For a given trait value \code{t}, \code{filt(t)} represents the probability that an individual with trait \code{t} enters the local community. If \code{filt = NULL},  \code{coalesc()} provides a neutral community.
}
 \item{add}{
indicates if additional variables must be passed to \code{filt}. It can be, for instance, environmental data conditioning the trait-based filtering in the community. Default is FALSE.
}
 \item{var.add}{
additional variables to be passed to \code{filt} when \code{add = T}.
}
 \item{pool}{
the regional pool of species providing immigrants to the local community. It should include the label of individual on first column, and of its species on second column. If \code{pool = NULL}, the pool is simulated as a metacommunity at speciation-drift equilibrium, based on parameter \code{theta}. The provided \code{pool} can contain trait values for each individuals in a third column.
}
 \item{traits}{
a matrix or data.frame including one or several traits on columns. A unique trait value is assigned to each species in the regional pool. Species names of \code{pool} must be included in \code{traits}.
If \code{traits = NULL}  and trait information is absent from \code{pool}, a random trait value is given to species of the regional pool, from a uniform distribution between 0 and 1. 
}
 \item{Jpool}{
if \code{pool = NULL}, it is the number of individuals to be simulated in the regional pool.
}
 \item{verbose}{
if \code{verbose = TRUE}, functions returns a lot of outputs about parameters, species pool and environmental filter.
}
}
\details{
Coalescent-based simulation of a community of size \code{J}. This generic function can simulate a neutral community (if \code{filt = NULL}) or a community undergoing both neutral and niche-based dynamics. In the latter case, \code{filt(t)} represents the relative ability of immigrants with trait values \code{t} in the regional pool to enter the community.
}
\value{
  \item{com}{
a data.frame of simulated individuals, with the label of ancestor individual in the regional pool on first column (as in first column of pool), species label on second column (as in second column of pool), and species trait (as in third column of pool).
\strong{Not provided if \code{m = 1} and \code{filt = NULL}: in this case the function provides a sample of the regional pool.}}
  \item{pool}{
a data.frame of the individuals of the regional source pool, with the label of ancestor individual in the regional pool on first column (as in first column of input \code{pool}), species label on second column (as in second column of input \code{pool}), and species trait (as in third column of input \code{pool}).}
  \item{call}{
    the call function.
  }
}
\references{
Hurtt, G. C. and S. W. Pacala (1995). "The consequences of recruitment limitation: reconciling chance, history and competitive differences between plants." Journal of Theoretical Biology 176(1): 1-12.

Gravel, D., C. D. Canham, M. Beaudet and C. Messier (2006). "Reconciling niche and neutrality: the continuum hypothesis." Ecology Letters 9(4): 399-409.

Munoz, F., P. Couteron, B. R. Ramesh and R. S. Etienne (2007). "Estimating parameters of neutral communities: from one Single Large to Several Small samples." Ecology 88(10): 2482-2488.

Munoz, F., B. R. Ramesh and P. Couteron (2014). "How do habitat filtering and niche conservatism affect community composition at different taxonomic resolutions?" Ecology 95(8): 2179-2191. 
}
\author{
F. Munoz
}

\examples{
# Simulation of a neutral community including 100 individuals
J <- 500; theta <- 50; m <- 0.1
comm1 <- coalesc(J, m, theta) 
# Regional and local trait distributions
plot_comm(comm1)
 
# Define a regional pool of species with equal abundances
pool <- cbind(1:10000, rep(1:500, 20), rep(NA, 10000))
# Uniform distribution of trait values
t.sp <- runif(500)
# No intraspecific variation
pool[,3] <- t.sp[pool[,2]]
# Generate a neutral community drawn from the pool
comm2<- coalesc(J, m, pool = pool) 
plot_comm(comm2)

# Directional habitat filtering toward t = 0
comm3 <- coalesc(J, m, filt = function(x) 1 - x, pool = pool)
# Regional and local trait distributions
plot_comm(comm3)

# Function for environmental filtering
sigma <- 0.1
filt_gaussian <- function(t, x) exp(-(x - t)^2/(2*sigma^2))

# Stabilizing habitat filtering around t = 0.1
comm4a <- coalesc(J, m, filt = function(x) filt_gaussian(0.1, x), pool = pool)
plot_comm(comm4a)
# Stabilizing habitat filtering around t = 0.5
comm4b <- coalesc(J, m, theta, filt = function(x) filt_gaussian(0.5, x),
                  pool = pool)
plot_comm(comm4b)
# Stabilizing habitat filtering around t = 0.9
comm4c <- coalesc(J, m, theta, filt = function(x) filt_gaussian(0.9, x),
                  pool = pool)
plot_comm(comm4c)

# Mean trait values in communities reflect the influence of habitat filtering
mean(comm4a$com[, 3])
mean(comm4b$com[, 3])
mean(comm4c$com[, 3])

# Disruptive habitat filtering around t = 0.5
comm5 <- coalesc(J, m, filt = function(x) abs(0.5 - x), pool = pool)
plot_comm(comm5)

# Multi-modal habitat filtering 
t.sp <- rnorm(500)
pool[, 3] <- t.sp[pool[,2]]
comm6 <- coalesc(J, m, filt = function(x) sin(3*x) + 1, pool = pool)
plot_comm(comm6)

}

\keyword{coalescent}
\keyword{neutral dynamics}
\keyword{niche-based dynamics}
