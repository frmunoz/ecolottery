\name{tcor}
\alias{tcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generates Correlated Traits
}
\description{
Create two random vectors of traits correlated between each other or a vector
of traits correlated to an existing one. The strength of the correlation is
defined by the parameter \code{rho}.
}
\usage{
tcor(n, rho = 0.5, mar.fun = rnorm, x = NULL, ...)
}
\arguments{
 \item{n}{
 a numeric parameter describing the number of value to generated.
}
 \item{rho}{
 a numeric parameter defining the strength of correlation between the two
 traits (Default = 0.5).
 It has to belong strictly to the [0, 1] interval.
}
 \item{x}{
 a vector which describe an existing trait vector. Default to NULL. 
}
 \item{mar.fun}{
 a function giving the random generation of the trait distribution, uses
 \code{rnorm} by default.
 }
 \item{...}{
  other arguments needed by \code{mar.fun()} function.
 }
}
\details{
  rho parameter is set to 0.5 by default. \code{x = NULL} by default.
  Code adapted from: \url{http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable}
}
\value{
 Return a data.frame with two numeric columns, each column defining a trait.
}
\author{
 P. Denelle
}

\examples{
# With no predefined trait
traits <- tcor(n = 10000, rho = 0.8)
plot(traits[, 1], traits[, 2])
cor(traits[, 1], traits[, 2])

# With existing trait
existing_trait <- rnorm(10000, 10, 1)
traits <- tcor(n = 10000, rho = 0.8, x = existing_trait)
plot(traits[, 1], traits[, 2])
cor(traits[, 1], traits[, 2])
}

\keyword{trait distribution}
\keyword{trait correlation}
