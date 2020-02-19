\name{forward}
\alias{forward}
\alias{gauss_limit}
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
forward(initial, m = 1, theta = NULL, d = 1, gens = 150, keep = FALSE,
        pool = NULL, traits = NULL, filt = NULL, filt.vect = F, limit.sim = NULL,
        limit.intra = F, par.limit = 0.1, coeff.lim.sim = 1, 
        type.filt = "immig", type.limit = "death", add = F, var.add = NULL, 
        prob.death = NULL, method.dist = "euclidean", checks = T, plot_gens = FALSE)
gauss_limit(dist, par)
get_number_of_gens(given_size, pool, traits = NULL, nbrep = 5, m = 1, theta = NULL,
                   d = 1, gens = NULL, filt = NULL,  limit.sim = NULL,
                   par.limit = 0.1, coeff.lim.sim = 1, type.filt = "immig", 
                   type.limit = "death", add = F, var.add = NULL, 
                   prob.death = NULL, method.dist = "euclidean", 
                   plot_gens = FALSE)
pick(com, d = 1, m = 1, theta = NULL, pool = NULL, prob.death = NULL,
     filt = NULL, filt.vect = F, limit.sim = NULL,  limit.intra = F, par.limit = 0.1, 
     coeff.lim.sim = 1, type.filt = "immig", type.limit = "death", add = add, 
     var.add = var.add, new.index = new.index, method.dist = "euclidean")
pick.mutate(com, d = 1, mu = 0, new.index = 0)
pick.immigrate(com, d = 1, m = 1, pool, prob.death = NULL, 
              filt = NULL, filt.vect = F, limit.sim = NULL,  limit.intra = F, 
              par.limit = 0.1, coeff.lim.sim = 1, 
              type.filt = "immig", type.limit = "death", add = F, var.add = NULL,
              method.dist = "euclidean")
}

\arguments{
 \item{com, initial}{
starting community. It is in principle a three (or more) column matrix or
data.frame including individual ID, species names and trait values. For strictly
neutral dynamics, it can be a vector of individual species names.
}
 \item{m}{
migration rate (if \code{m = 1} the community is a subsample of the regional pool).
}
 \item{theta}{
parameter of neutral dynamics in the regional pool (used only if \code{pool=NULL}), it is the \dQuote{fundamental biodiversity number} (\eqn{\theta}).
}
 \item{mu}{
mutation rate derived from \eqn{\theta}, such as \eqn{\mu = \theta/(2*J_m), where \eqn{J_m} is the size of the metacommunity.}.
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
a matrix or data.frame including one or several traits on columns. A unique trait value is assigned to each species in the regional pool. Species names of \code{pool} must be included in rownames of \code{traits}.
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
if non null, limiting similarity will be simulated, based on species trait
distances (computed with the method given by \code{method.dist}) and a Gaussian
overlapping function. 
}
 \item{limit.intra}{
should limiting similarity play on individual dynamics within species? Default is \code{FALSE}. 
}
 \item{par.limit}{
a vector of additional parameters to be passed to the limiting similarity function defined in \code{limit.sim}.
If \code{NULL} but the \code{limit.sim} requires additional parameter values,
an error message will be issues.
In the case of default Gaussian limiting similarity function, the additional 
parameter value sets the standard deviation of the Gaussian function.
}
 \item{coeff.lim.sim}{
adjust the intensity of limiting similarity. 
}
 \item{type.filt}{
indicates how habitat filtering plays in community assembly. 
If \code{type.filt="death"}, it plays on mortality events.
If \code{type.filt="immig"} (default), it plays at the establishment of immigrants,
If \code{type.filt="loc.recr"}, it plays at the establishment of local offspring.
\code{type.filt} can be any combination of these three values.
}
\item{type.limit}{
indicates how limiting similarity plays in community assembly. 
If \code{type.limit="death"} (default), it plays on mortality events.
If \code{type.limit="immig"}, it plays at the establishment of immigrants,
If \code{type.limit="loc.recr"}, it plays at the establishment of local offspring.
\code{type.limit} can be any combination of these three values.
}
 \item{filt}{
the function used to represent habitat filtering. For a given trait value
\code{t}, \code{filt(t)} represents the probability that an individual with
trait \code{t} enters the local community.
}
 \item{filt.vect}{
indicates whether the filtering function can be vectorized. It means that the function can take as input a vector of trait values and provide a vector of the corresponding weights.
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
provides the method to compute trait distances between individuals (syntax of
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
 \item{checks}{
should initial checks that the inputs are correct be performed?
}
 \item{dist}{
vector of trait distances that is used for calculation of limiting similarity.
}
 \item{par}{
addditional parameter values to be passed to \code{limit.sim} function.
}
}
\details{
The model of community assembly is a zero-sum game, so that the number of individuals 
of the community is fixed to the number of individuals in initial community.

Two types of niche-based processes, habitat filtering and limiting similarity
can affect (i) immigration, (ii) mortality, and (iii) recruitment of local offspring.
The user can define on each of the three components the processes play, by defining 
\code{type.filt} and  \code{type.limit}. These parameters can combine one of several
values in "immig", "death" and "loc.recr" depending on either the process plays on 
immigration, (ii) mortality, and (iii) recruitment of local offspring, respectively.
For instance, defining \code{type.filt = c("immig", "death")} would mean that habitat 
filtering plays on both immigration and on death events. By default, \code{type.filt = "immig"}
and \code{type.simil = "death"}, which corresponds to a classical conception of 
the hierarchical influence of habitat filtering and limiting similarity on immigration success
and local survival.

A environmental filtering function can be defined with \code{filt}. It should take two arguments, the first is a trait value of a candidate immigrant, the second is a vector including the parameter values of the filtering function. See examples below for further information on usage.
An important point is that in many cases the function might be vectorized, that is, it can provide a vector of filtering probabilities for a vector of trait values given in first argument. In this case the user should set \code{filt.vect = T}, which will significantly accelerate simulations. 

\code{gauss_limit()} is the default function used to compute limiting similarity. In this case, the relative performance of an individual decreases depending on the sum of exponential distances to other individuals in the community (McArthur and Levins 1967), i.e.,

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

A \href{https://ecolottery.shinyapps.io/forward/}{Shiny app} is available to visualize simulated trait distributions for chosen parameter values in the model.
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
For neutral dynamics, see S. P. Hubbell 2001. "The Unified Neutral Theory of
Biodiversity". Princeton University Press. 

For the default model of limiting similarity,
MacArthur, R., Levins, R., 1967. The limiting similarity, convergence, and divergence of coexisting species. The American Naturalist 101, 377–385.
}
\author{
F. Munoz, derived from the \code{untb} function of R. Hankin.
}

\examples{
\dontrun{
# Initial community composed of 10 species each including 10 individuals
initial1 <- rep(as.character(1:10), each = 10)

# Simulation of speciation and drift dynamics over 1000 time steps
final1 <- forward(initial = initial1, theta = 50, gens = 1000)
# The final community includes new species (by default names begins with "new")
final1$com$sp # includes new species generated by speciation events

# A regional pool including 100 species each including 10 individuals
pool <- rep(paste("pool.sp",as.character(1:100), sep=""), each = 10)

# Simulation of migration and drift dynamics over 1000 time steps
final2 <- forward(initial = initial1, m = 0.1, gens = 1000, pool = pool)
# The final community includes species that have immigrated from the pool
final2$com$sp # includes new species that immigrated from the pool

# Initial community composed of 10 species each including 10 individuals, 
# with trait information for niche-based dynamics
initial2 <- data.frame(ind = 1:100, sp = rep(as.character(1:10), each = 10), 
                      trait = runif(100), stringsAsFactors = F)

# Simulation of stabilizing hab. filtering around t = 0.5, over 2000 time steps
sigm <- 0.1
filt_gaussian <- function(t,x) exp(-(x - t)^2/(2*sigm^2))
# Filtering only plays only on immigration (default type.filt = "immig")
final3 <- forward(initial = initial2, m = 0.1, gens = 2000, pool = pool, 
                 filt = function(x) filt_gaussian(0.5,x))
plot_comm(final3) # trait distribution in final community
# Filtering also plays on local mortality and recruitment
final4 <- forward(initial = initial2, m = 0.1, gens = 2000, pool = pool, 
                 filt = function(x) filt_gaussian(0.5,x), 
                 type.filt = c("immig", "death", "loc.recr"))
plot_comm(final4) # Stronger convergence around 0.5 

# Simulation of stabilizing hab. filtering with two traits
pool <- data.frame(ind=1:10000, sp=rep(1:500, 20), tra1=rep(NA, 10000), tra2=rep(NA, 10000))
# Trait values in the pool
t1.sp <- runif(500)
t2.sp <- runif(500)
pool[, 3] <- t1.sp[pool[,2]]
pool[, 4] <- t2.sp[pool[,2]]
initial3 <- pool[sample(1:10000, 100),]
# Parameters of filtering for first trait
topt1 <- 0.25
sigm1 <- 0.1
# Parameters of filtering for second trait
topt2 <- 0.75
sigm2 <- 0.2
filt_gaussian <- function(x) exp(-(x[1] - topt1)^2/(2*sigm1^2)) * exp(-(x[2] - topt2)^2/(2*sigm2^2))
system.time(final5 <- forward(initial = initial3, m = 0.1, gens = 2000, pool = pool, 
                 filt = filt_gaussian))
# Here the filtering function can be vectorized
system.time(final5 <- forward(initial = initial3, m = 0.1, gens = 2000, pool = pool, 
                 filt = filt_gaussian, filt.vect = T))
plot_comm(final5, seltrait = 1) # trait distribution in final community (first trait)
plot_comm(final5, seltrait = 2) # trait distribution in final community (second trait)

# Simulation of limiting similarity, over 2000 time steps
final6 <- forward(initial = initial2, m = 0.1, gens = 2000, pool = pool, 
                 limit.sim = TRUE)
plot_comm(final6)
# Check temporal changes
plot(final6$sp_t, xlab = "Time step", ylab = "Community richness")
# Index of limiting similarity over time
plot(final6$dist.t, xlab = "Time step", ylab = "Limiting similarity")

# Higher migration rate, 5000 time steps
# Limiting similarity only plays on mortality
final7 <- forward(initial = initial2, m = 0.7, gens = 5000, pool = pool, 
                 limit.sim = TRUE)
plot_comm(final7) # should be closer to regional distribution
# Limiting similarity plays on immigration, mortality and local recruitment
final8 <- forward(initial = initial2, m = 0.7, gens = 5000, pool = pool, 
                 limit.sim = TRUE, type.limit = c("immig", "death", "loc.recr"))
plot_comm(final8)
# Check temporal changes
plot(final8$sp_t, xlab = "Time step", ylab = "Community richness")
points(1:5000, final7$sp_t, col="blue")
# Index of limiting similarity over time
plot(final8$dist.t, xlab = "Time step", ylab = "Limiting similarity")
points(1:5000, final7$dist.t, col="blue")
}
}

\keyword{neutral dynamics}
\keyword{niche-based dynamics}
