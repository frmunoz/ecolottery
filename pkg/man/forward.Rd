\name{forward}
\alias{forward}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of neutral and niche-based community dynamics forward in time
}
\description{
Simulates niche-based (habitat filtering and/or limiting similarity) and neutral community dynamics from a given initial composition, over a given number of generations. 
}
\usage{
forward(initial, prob = 0, D = 1, gens = 150, keep = FALSE, pool = NULL, limit.sim = F, coeff.lim.sim = 2, sigma=0.1, filt = NULL, prob.death = NULL, method.dist = "euclidean")
select (a, D = 1, prob = 0, pool = NULL, prob.death = prob.death, limit.sim = NULL, coeff.lim.sim = 2, sigma = 0.1, filt = NULL, new.index = new.index)
select.mutate(a, D = 1, prob.of.mutate = 0, new.index = 0)
select.immigrate (a, D = 1, prob.of.immigrate = 0,  pool, prob.death= NULL, limit.sim = NULL, coeff.lim.sim = 2, sigma = 0.1, filt = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a, initial}{
Starting community. It is in principle a three-column matrix or data frame including individual ID, species names and trait values. For strictly neutral dynamics, can be a vector of individual species names.
}
  \item{prob, prob.of.immigrate, prob.of.mutate}{
Probability of an individual establishing in the community not being a descendent of an existing individual. If descendent from a new ancestor, can be either through immigration.
}
  \item{D}{
Number of individuals that die in each time-step.
}
  \item{gens}{
Number of generations to simulate.
}
  \item{keep}{
Boolean with default FALSE meaning to return community composition at the end of the simulation, and TRUE meaning to return a list including community composition at successive time steps (see Value section).
}
  \item{pool}{
The regional pool of species providing immigrants to the local community. It is in principle a three-column matrix or data frame including individual ID, species names and trait values. If trait information is missing, a random trait value is given to individuals, from a uniform distribution between 0 and 1. If NULL, the pool is simulated as a metacommunity at speciation-drift equilibrium, based on prob for speciation rate.
}
  \item{limit.sim}{
If T, limiting similarity will be simulated, based on species trait distances (computed with method given in method.dist). 
}
  \item{coeff.lim.sim}{
Adjust the intensity of limiting similarity. 
}
  \item{sigma}{
Adjust the variance of the overlap function used to calculate limiting similarity.
}
 \item{filt}{
The function used to represent habitat filtering. For a given trait value, t, filt(t) represents the probability that an individual with trait t enters the local community.
}
  \item{prob.death}{
Provides a baseline probability of death that is homogeneous across species. It is used in niche-based dynamics to represent the balance of baseline and niche-dependent mortality.
}
  \item{method.dist}{
Provide the method to compute trait distances between individuals (syntax of function dist).
}
  \item{new.index}{
Used to give a new species name when speciation occurs.
}
}
\details{
When niche-based dynamics are simulated, the niche-based constraints influence both immigration and mortality.

It is a zero-sum game, so that the number of individuals of the community is fixed to the number of individuals in initial community.

Functions select.immigrate() and select.mutate() are used to simulate immigration and speciation events within a time step. They are embedded in forward and are not really intended for the end user. 
}
\value{
\item{a}{
A matrix of individual in the final community, with individual ID in first column, species name in second column, and trait value in the third column.
}
\item{pool}{
The regional pool of species used in simulation. Includes label of individual on first column, its species on second column, and a trait attribute on the third column.
}
\item{aa}{
If keep = T, a list of community composition for each time step (a matrix as in a).
}
\item{limit.sim.t}{
If limit.sim = T, the average value of the limiting similarity function over time.
}
\item{new.index}{
For  select.mutate(), return the new index to be used for species name at a next speciation event.
}
}
\references{
For neutral dynamics, S. P. Hubbell 2001. “The Unified Neutral Theory of Biodiversity”. Princeton University Press. 
}
\author{
François Munoz, derived from the untb function of Robin K. S. Hankin
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Initial community composed of 10 species each including 10 individuals
initial <- sort(rep(as.character(1:10),10))

# Simulation of speciation and drift dynamics over 100 time steps
final <- forward(initial=initial, prob=0.1, gens=100)
# The final community includes new species (by default their name begins with "new")
final$a[,2] # includes new species generated by speciation events

# A regional pool including 100 species each including 10 individuals
pool <- sort(rep(as.character(1:100),10))

# Simulation of migration and drift dynamics over 100 time steps
final <- forward(initial=initial, prob=0.1, gens=100, pool=pool)
# The final community includes species that have immigrated from the pool
final$a[,2] # includes new species that immigrated from the pool

# Initial community composed of 10 species each including 10 individuals, with trait information for niche-based dynamics
initial <- data.frame(sp=sort(rep(as.character(1:10),10)),trait=runif(100))

# Simulation of stabilizing habitat filtering around t=0.5, over 100 time steps
final <- forward(initial=initial, prob=0.1, gens=100, pool=pool, filt=function(x) 0.5-abs(0.5-x))
hist(final$a[,3]) # histogram of trait values in final community
# With higher immigration
final <- forward(initial=initial, prob=0.8, gens=100, pool=pool, filt=function(x) 0.5-abs(0.5-x))
hist(final$a[,3]) # should be closer to 0.5

# Simulation of limiting similarity, over 100 time steps
final <- forward(initial=initial, prob=0.1, gens=100, pool=pool, limit.sim=T)
hist(final$a[,3])
# Stronger limiting similarity
final <- forward(initial=initial, prob=0.1, gens=100, pool=pool, limit.sim=T, coeff.lim.sim = 20)
hist(final$a[,3]) # the distribution will be more even

# Variation of community richness with time
final <- forward(initial=initial, prob=0.1, gens=1000, pool=pool, limit.sim=T, keep=T)
plot(unlist(lapply(final$aa,function(x) length(unique(x[,2])))),xlab="Time step",ylab="Community richness") # use it to check stationarity

# Index of limiting similarity over time
plot(final$ limit.sim.t,xlab="Time step",ylab="Limiting similarity")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ neutral dynamics }
\keyword{ niche-based dynamics }
