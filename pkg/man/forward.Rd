\name{forward}
\alias{forward}
\alias{get_number_of_gens}
\alias{pick}
\alias{pick.mutate}
\alias{pick.immigrate}
\title{
Simulation of neutral and niche-based community dynamics forward in time
}
\description{
Simulates niche-based (habitat filtering and/or limiting similarity) and neutral
community dynamics from a given initial composition, over a given number of
generations. 
}

\usage{
forward(initial, prob = 0, d = 1, gens = 150, keep = FALSE,
        pool = NULL, traits = NULL, limit.sim = FALSE, coeff.lim.sim = 1,
        type.assemb = "death", sigm = 0.1, filt = NULL, add = F, var.add = NULL, 
        prob.death = NULL, method.dist = "euclidean", plot_gens = FALSE)
get_number_of_gens(given_size, pool, traits = NULL, nbrep = 5, prob = 1, d = 1, 
                   gens = NULL, limit.sim = FALSE, 
                   coeff.lim.sim = 1, type.assemb = "death", sigm = 0.1, filt = NULL, 
                   add = F, var.add = NULL, prob.death = NULL, 
                   method.dist = "euclidean", plot_gens = FALSE)
pick(com, d = 1, prob = 0, pool = NULL, prob.death = prob.death,
     limit.sim = NULL, coeff.lim.sim = 1, type.assemb = "death", sigm = 0.1, filt = NULL, 
     new.index = new.index, method.dist = "euclidean")
pick.mutate(com, d = 1, prob.of.mutate = 0, new.index = 0)
pick.immigrate(com, d = 1, prob.of.immigrate = 0, pool, prob.death = NULL, 
              limit.sim = NULL, coeff.lim.sim = 1, type.assemb = "death", sigm = 0.1, 
              filt = NULL, method.dist = "euclidean")
}

\arguments{
 \item{com, initial}{
starting community. It is in principle a three (or more) column matrix or
data.frame including individual ID, species names and trait values. For strictly
neutral dynamics, it can be a vector of individual species names.
}
 \item{prob, prob.of.immigrate, prob.of.mutate}{
probability of an individual establishing in the community not being a
descendant of an existing individual. If descendant from a new ancestor, can be
either through immigration (in \code{pick.immigrate()}) or through mutation
(in \code{pick.mutate()}).
}
 \item{d}{
number of individuals that die in each time-step.
}
 \item{gens}{
number of generations to simulate.
}
 \item{keep}{
boolean value. If \code{FALSE} (default) the function output only the community
composition at the end of the simulation. If \code{TRUE} the function output a
list of community composition at successive time steps (see Value section).
}
 \item{pool}{
the regional pool of species providing immigrants to the local community. It is
in principle a three-column matrix or data frame including individual ID,
species names and trait values. If trait information is missing, a random trait
value is given to individuals, from a uniform distribution between 0 and 1.
If \code{NULL}, the pool is simulated as a metacommunity at speciation-drift
equilibrium, based on \code{prob} for speciation rate.
}
 \item{traits}{
a matrix or data.frame including one or several traits on columns. A unique trait value is assigned to each species in the regional pool. Species names of \code{pool} must be included in \code{traits}.
If \code{traits = NULL} and trait information is absent from \code{pool}, a random trait value is given to species of the regional pool, from a uniform distribution between 0 and 1. 
}
 \item{given_size}{
size of the community you want to have an estimate of the number of generations
needed to reach stationarity in species richness.
}
\item{nbrep}{
number of replicates from which you want to estimate the number of generations
needed to reach stationarity in species richness.
}
 \item{limit.sim}{
if \code{TRUE,} limiting similarity will be simulated, based on species trait
distances (computed with the method given by \code{method.dist}) and a Gaussian
overlapping function. 
}
 \item{coeff.lim.sim}{
adjust the intensity of limiting similarity. 
}
 \item{type.assemb}{
indicate how limiting similarity and/or habitat filtering play in community assembly. 
If \code{type.assemb="death"} (default), it plays on mortality events.
If \code{type.assemb="estab"}, it plays at the establishment of immigrants
and local offspring.
If \code{type.assemb="both"}, it plays at both mortality and establishment.
If \code{type.assemb="immig"}, it plays only at the establishment of immigrants.
}
 \item{sigm}{
adjust the variance of the overlap function used to calculate limiting
similarity.
}
 \item{filt}{
the function used to represent habitat filtering. For a given trait value
\code{t}, \code{filt(t)} represents the probability that an individual with
trait \code{t} enters the local community.
}
 \item{add}{
indicates if additional variables must be passed to \code{filt}. It can be, for instance, environmental data conditioning the trait-based filtering in the community. Default is FALSE.
}
 \item{var.add}{
additional variables to be passed to \code{filt} when \code{add = T}.
}
 \item{prob.death}{
provides a baseline probability of death that is homogeneous across species. It
is used in niche-based dynamics to represent the balance of baseline and
niche-dependent mortality.
}
 \item{method.dist}{
provide the method to compute trait distances between individuals (syntax of
function \code{\link[stats]{dist}}, can be in the list \code{c("euclidean",
"maximum", "manhattan", 
"canberra", "binary", "minkowski")}).
}
 \item{new.index}{
prefix used to give a new species name when speciation occurs.
}
 \item{plot_gens}{
plot the number of unique individuals and species over generations.
}
}
\details{
It is a zero-sum game, so that the number of individuals of the community is
fixed to the number of individuals in initial community.

When niche-based dynamics are simulated, the niche-based constraints influence
both immigration and mortality.

Function \code{get_number_of_gen()} allows determining the number of generations
needed to reach stationary richness for given parameterization of
\code{forward()}. The target number of generation is based on assessing the
change point in species richness change over time for replicate simulated
communities with random initial composition. A conservative measure is proposed
as the maximum time to reach stationary richness over the replicate simulated
communities.

Functions \code{pick.immigrate()} and \code{pick.mutate()} are used to simulate
immigration and speciation events within a time step. They are embedded in
forward and are not really intended for the end user. 
}
\value{
\item{com}{
if \code{keep = FALSE}, a data.frame of simulated individuals, with the label
of ancestor individual in the regional pool on first column (as in the first
column of the pool), species label on second column (as in the second column of
the pool), and species trait (as in the third column of the pool).
}
\item{pool}{
a data.frame of the individuals of the regional source pool, with the label of
ancestor individual in the regional pool on first column (as in first column of
input \code{pool}), species label on second column (as in second column of input
\code{pool}), and species trait (as in third column of input \code{pool}).
}
\item{sp_t}{
a vector of species richness at each time step.
}
\item{com_t}{
if \code{keep = TRUE}, a list of community composition for each time step 
(a data.frame as in \code{com}).
}
\item{dist.t}{
if \code{limit.sim = TRUE}, the average value of the limiting similarity
function over time.
}
\item{new.index}{
for \code{pick.mutate()}, return the new index to be used for species name at a
next speciation event.
}
\item{call}{
the call function.
}
}
\references{
For neutral dynamics, S. P. Hubbell 2001. "The Unified Neutral Theory of
Biodiversity". Princeton University Press. 
}
\author{
F. Munoz, derived from the \code{untb} function of R. Hankin.
}

\examples{
\dontrun{
# Initial community composed of 10 species each including 10 individuals
initial1 <- rep(as.character(1:10), each = 10)

# Simulation of speciation and drift dynamics over 1000 time steps
final1 <- forward(initial = initial1, prob = 0.1, gens = 1000)
# The final community includes new species (by default names begins with "new")
final1$com$sp # includes new species generated by speciation events

# A regional pool including 100 species each including 10 individuals
pool <- rep(paste("pool.sp",as.character(1:100), sep=""), each = 10)

# Simulation of migration and drift dynamics over 1000 time steps
final2 <- forward(initial = initial1, prob = 0.1, gens = 1000, pool = pool)
# The final community includes species that have immigrated from the pool
final2$com$sp # includes new species that immigrated from the pool

# Initial community composed of 10 species each including 10 individuals, 
# with trait information for niche-based dynamics
initial2 <- data.frame(sp = rep(as.character(1:10), each = 10), 
                      trait = runif(100))

# Simulation of stabilizing hab. filtering around t = 0.5, over 2000 time steps
sigm <- 0.1
filt_gaussian <- function(t,x) exp(-(x - t)^2/(2*sigm^2))
# Filtering only plays on mortality (default type.assemb = "death")
final3 <- forward(initial = initial2, prob = 0.1, gens = 2000, pool = pool, 
                 filt = function(x) filt_gaussian(0.5,x))
plot_comm(final3) # trait distribution in final community
# Filtering also plays on establishment (default type.assemb = "both")
final4 <- forward(initial = initial2, prob = 0.1, gens = 2000, pool = pool, 
                 filt = function(x) filt_gaussian(0.5,x), type.assemb = "both")
plot_comm(final4) # Stronger convergence around 0.5 

# Simulation of limiting similarity, over 1000 time steps
final5 <- forward(initial = initial2, prob = 0.1, gens = 1000, pool = pool, 
                 limit.sim = TRUE)
plot_comm(final5)
# Check temporal changes
plot(final5$sp_t, xlab = "Time step", ylab = "Community richness")
# Index of limiting similarity over time
plot(final5$dist.t, xlab = "Time step", ylab = "Limiting similarity")

# Higher migration rate
# Without limiting similarity on establishment
final6 <- forward(initial = initial2, prob = 0.8, gens = 1000, pool = pool, 
                 limit.sim = TRUE)
plot_comm(final6) # should be closer to regional distribution
# With limiting similarity on establishment
final7 <- forward(initial = initial2, prob = 0.8, gens = 1000, pool = pool, 
                 limit.sim = TRUE, type.assemb = "both")
plot_comm(final7)
# Check temporal changes
plot(final7$sp_t, xlab = "Time step", ylab = "Community richness")
points(1:1000, final6$sp_t, col="blue")
# Index of limiting similarity over time
plot(final7$dist.t, xlab = "Time step", ylab = "Limiting similarity")
points(1:1000, final6$dist.t, col="blue")
}
}

\keyword{neutral dynamics}
\keyword{niche-based dynamics}
