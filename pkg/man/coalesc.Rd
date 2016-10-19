\name{coalesc}
\alias{coalesc}
\alias{coalesc_simp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coalescent-based simulation of ecological communities undergoing both neutral and niche-base dynamics
}
\description{
Simulates the composition of a community based on immigration from a regional pool, habitat filtering depending on local environment and species traits, and local birth-death stochastic dynamics. 
}
\usage{
coalesc(J, m = 1, theta = NULL, filt = NULL, pool = NULL, traits = NULL, Jpool = 50 * J)
coalesc_simp(J, m = 1)
}
\arguments{
 \item{J}{
Number of individuals in the local community.
}
 \item{m}{
Migration rate (if m=1 the community is a subsample of the regional pool).
}
 \item{theta}{
Parameter of neutral dynamics in the regional pool (used only if pool=NULL).
}
 \item{filt}{
A function representing the effect of local habitat filtering. For a given trait value, t, filt(t) represents the probability that an individual with trait t enters the local community. If filt=NULL, the lottery function provides a neutral community.
}
 \item{pool}{
The regional pool of species providing immigrants to the local community. It should include the label of individual on first column, and of its species on second column. If NULL, the pool is simulated as a metacommunity at speciation-drift equilibrium, based on parameter theta.
}
 \item{traits}{
A matrix or data.frame including one or several traits on columns. In the present version of the function, a unique trait value is assigned to each species un the regional pool. If NULL, a random trait value is given to species of the regional pool, from a uniform distribution between 0 and 1.
}
 \item{Jpool}{
If pool = NULL, represents the number of individuals to be simulated in the regional pool.
}
}
\details{
Coalescent-based simulation of a community of size J. This generic function can simulate a neutral community (if filt=NULL) or a community undergoing both neutral and niche-based dynamics. In the later case, filt(t) represents the relative ability of immigrants with trait values t in the regional pool to enter the community.
.coalesc is a simplified version of coalesc to increase speed in intensive calculations.
}
\value{
\item{com}{
A matrix of simulated individuals, with the label of ancestor individual in the regional pool on first column (as in first column of pool), species label on second column (as in second column of pool), and species trait (as in third column of pool).
Not provided if m=1 and filt=NULL: in this case the function provides a sample of the regional pool.
}
\item{pool}{
A matrix of the individuals of the regional source pool, with the label of ancestor individual in the regional pool on first column (as in first column of pool), species label on second column (as in second column of pool), and species trait (as in third column of pool).
}
}
\references{
Hurtt, G. C. and S. W. Pacala (1995). "The consequences of recruitment limitation: reconciling chance, history and competitive differences between plants." Journal of Theoretical Biology 176(1): 1-12.

Gravel, D., C. D. Canham, M. Beaudet and C. Messier (2006). "Reconciling niche and neutrality: the continuum hypothesis." Ecology Letters 9(4): 399-409.

Munoz, F., P. Couteron, B. R. Ramesh and R. S. Etienne (2007). "Estimating parameters of neutral communities: from one Single Large to Several Small samples." Ecology 88(10): 2482-2488.

Munoz, F., B. R. Ramesh and P. Couteron (2014). "How do habitat filtering and niche conservatism affect community composition at different taxonomic resolutions?" Ecology 95(8): 2179-2191. 
}
\author{
Fran\c{c}ois Munoz
}

\examples{
# Simulation of a neutral community including 100 individuals
J <- 500; theta <- 50; m <- 0.1;
comm1 <- coalesc(J, m, theta) 
# Regional and local trait distributions
plot_comm(comm1)
 
# Define a regional pool of species with equal abundances
pool <- cbind(1:10000, rep(1:500,20), rep(NA,10000))
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
plot_comm(comm2)

# Function for environmental filtering
sigma <- 0.1;
filt_gaussian <- function(t,x) exp(-(x-t)^2/(2*sigma^2))

# Stabilizing habitat filtering around t = 0.1
comm4a <- coalesc(J, m, filt = function(x) filt_gaussian(0.1,x), pool = pool)
plot_comm(comm4a)
# Stabilizing habitat filtering around t = 0.5
comm4b <- coalesc(J, m, theta, filt = function(x) filt_gaussian(0.5,x), pool = pool)
plot_comm(comm4b)
# Stabilizing habitat filtering around t = 0.9
comm4c <- coalesc(J, m, theta, filt = function(x) filt_gaussian(0.9,x), pool = pool)
plot_comm(comm4c)
# Mean trait values in communities reflect the influence of habitat filtering
mean(comm4a$com[, 3])
mean(comm4b$com[, 3])
mean(comm4c$com[, 3])

# Disruptive habitat filtering around t = 0.5
comm5 <- coalesc(J, m, filt = function(x) abs(0.5 - x), pool = pool)
plot_comm(comm5)

}

\keyword{coalescent}
\keyword{neutral dynamics}
\keyword{niche-based dynamics}
