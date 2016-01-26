\name{coalesc}
\alias{coalesc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coalescent-based simulation of ecological communities undergoing both neutral and niche-base dynamics
}
\description{
Simulates the composition of a community based on immigration from a regional pool, habitat filtering depending on local environment and species traits, and local birth-death stochastic dynamics. 
}
\usage{
coalesc(J, theta, m = 1, filt = NULL, pool = NULL, Jpool = 50 * J)
}
\arguments{
 \item{J}{
Number of individuals in the local community.
}
 \item{theta}{
Parameter of neutral dynamics in the regional pool (used only if pool=NULL).
}
 \item{m}{
Migration rate (if m=1 the community is a subsample of the regional pool).
}
 \item{filt}{
A function representing the effect of local habitat filtering. For a given trait value, t, filt(t) represents the probability that an individual with trait t enters the local community. If filt=NULL, the lottery function provides a neutral community.
}
 \item{pool}{
The regional pool of species providing immigrants to the local community. If NULL, the pool is simulated as a metacommunity at speciation-drift equilibrium, based on parameter theta. In this case, a random trait value is given to individuals of the regional pool, from a uniform distribution between 0 and 1. Otherwise, input value of pool must include the label of individual on first column, of its species on second column, and a trait attribute on the third column.
}
 \item{Jpool}{
If pool = NULL, represents the number of individuals to be simulated in the regional pool.
}
}
\details{
Coalescent-based simulation of a community of size J. This generic function can simulate a neutral community (if filt=NULL) or a community undergoing both neutral and niche-based dynamics. In the later case, filt(t) represents the relative ability of immigrants with trait values t in the regional pool to enter the community.
}
\value{
\item{res}{
A matrix of simulated individuals, with the label of ancestor individual in the regional pool on first column (as in first column of pool), species label on second column (as in second column of pool), and species trait (as in third column of pool).
}
}
\references{
Hurtt, G. C. and S. W. Pacala (1995). "The consequences of recruitment limitation: reconciling chance, history and competitive differences between plants." Journal of Theoretical Biology 176(1): 1-12.

Gravel, D., C. D. Canham, M. Beaudet and C. Messier (2006). "Reconciling niche and neutrality: the continuum hypothesis." Ecology Letters 9(4): 399-409.

Munoz, F., P. Couteron, B. R. Ramesh and R. S. Etienne (2007). "Estimating parameters of neutral communities: from one Single Large to Several Small samples." Ecology 88(10): 2482-2488.

Munoz, F., B. R. Ramesh and P. Couteron (2014). "How do habitat filtering and niche conservatism affect community composition at different taxonomic resolutions?" Ecology 95(8): 2179-2191. 
}
\author{
François Munoz
}

\examples{
# Simulation of a neutral community including 100 individuals
J <- 100; theta <- 50; m <- 0.1;
comm1 <- coalesc(J, theta, m) 
hist(comm1[, 3])

# Directional habitat filtering toward t = 0
comm2 <- coalesc(J, theta, m, filt = function(x) 1 - x)
hist(comm2[, 3])

# Stabilizing habitat filtering around t = 0.5
comm3 <- coalesc(J, theta, m, filt = function(x) 0.5 - abs(0.5 - x))
hist(comm3[, 3])

# Disruptive habitat filtering around t = 0.5
comm4 <- coalesc(J, theta, m, filt = function(x) abs(0.5 - x))
hist(comm4[, 3])
}


\keyword{coalescent}
\keyword{neutral dynamics}
\keyword{niche-based dynamics}
