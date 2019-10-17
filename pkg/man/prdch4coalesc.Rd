\name{prdch4coalesc}
\alias{prdch4coalesc}

\title{
Posterior predictive checks using coalescent-based simulation of ecological communities
}
\description{
This function performs posterior predicitve checks assessing a model's prediction performance of either user defined summary statistics, functional diversity (moments of the trait distribution, Enquist et al., 2015) or species rank-abundance plots}

\usage{
prdch4coalesc <- function(com.obs, pool, filt, params, stats = "abund", f.stats = NULL, estim = NULL)
}

\arguments{
  \item{com.obs}{
a data.frame of observed individuals, the first column should contain individual labels, the second the names of the species and the third species trait values - should be the same as provided for the ABC analysis in \code{coalesc-abc}}
  
  \item{pool}{
a data.frame containing the regional pool of species providing immigrants to the local community. It should include the label of individual on first column, its species on second column and their trait values in the third column as in \code{com.obs} - shuld be the same as provided for the ABC analysis in \code{coalesc-abc}}

  \item{filt}{
the inferred environmental filtering function. If \code{filt = NULL} predicitive checks will be performed for a purely neutral model and \code{params} should only be one column containing the posterior distribution of the migration rate parameter \code{m}}

  \item{params}{
the posterior distribution resulting from the ABC analysis, should be a data.frame containing each parameter distribution as columns (posterior distribution of the migration rate parameter (\code{m}) should be the last column) }
  
  \item{stats}{
statistics used for the predictive checks:
    \itemize{
      \item if \code{stats = "custom"}, predictive checks using a user-defined function provided in \code{f.stats}
      
      \item if \code{multi = "moments"}, predictive checks using the 4 first moments of the trait distribution evaluated by comparing the observed moments to those of 100 simulated communities under the model of interest
      
      \item if \code{multi = "abund"}, returns the confidence interval for the rank-abundance curves simulated by the model as well as the observed species rank-abundance curve distinguishing between species whose relative abundance is either over or under-estimated by the model}}
  
  \item{f.stats}{
the user-defined function computing the summary statistics to be used in the predictive checks when \code{stats = "custom"}}

  \item{estim}{
an estimator of the posterior distribution - can be either \code{"mean"}, \code{"median"} or \code{"mode"} used for the simulations} }

\value{
 \item{pvalue}{
    probabilities that the observed user-defined statistics are greater than
the same statistics simulated by the model (if \code{stats="custom"}) } 
    
 \item{prd.moments}{
    probabilities that the first four moments of the observed trait distribution are greater than those simulated by the model (if \code{stats="moments"})}
    
 \item{underestim.sp}{
    the names of the under-represented species by the model evaluated by comparing observed relative abundances of species to those simulated by the model (if \code{stats="abund"})}
    
 \item{overestim.sp}{
    the names of over-represented species by the model evaluated by comparing observed relative abundances of species to those simulated by the model (if \code{stats="abund"})} 

Also, if \code{stats="abund"}, rank-abundance curves are displayed distinguishing between species whose relative abundance is either over or under-estimated by the model as well as the confidence interval of the rank-abundance curves derived from the simulated communities }


\references{
  
  Enquist, Brian J., et al. "Scaling from traits to ecosystems: developing a general trait driver theory via integrating trait-based and metabolic scaling theories." Advances in ecological research. Vol. 52. Academic Press, 2015. 249-318.

Csillery, Katalin, et al. "Approximate Bayesian computation (ABC) in practice." Trends in ecology & evolution 25.7 (2010): 410-418.

}
\author{
E. Barthelemy
}

\examples{
pool <- data.frame(ind = 1:1000, sp = as.character(rep(1:50), each = 10),trait = runif(1000), stringsAsFactors = F)
com.obs <- pool[sample(nrow(pool), 50),]

f.stats <- function(com){
  tab <- table(com[,2])
  as.vector(t(sapply(0:3, function(x) hillR::hill_taxa(tab, q=x)))) }
  
filt <- function(x, par) exp(-(x - par[1])^2/(2*par[2]^2))

stats <- c("custom","abund", "moments")

params <- data.frame(par1 = runif(100,0,1), par2 = runif(100, .5, .9), par3 = runif(100,0,1))

checks <- prdch4coalesc(com.obs, pool, filt, params, stats, f.stats) }

\keyword{Coalescent}
\keyword{Approximate Bayesian Computation}
\keyword{Niche-based dynamics}
\keyword{Neutral dynamics}
\keyword{Predictive checks}
