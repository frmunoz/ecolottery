\name{plot.comm}
\alias{plot.comm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Regional vs. Local trait distrbutions of Regional vs. Local Abundances
}
\description{
  If type = "trait", the function provides density plots of the trait or abundance distributions in the regional pool and in a local community.
  If type = "abund", the function displays the relationship between regional and local species abundances.
  Default type is "trait".
  To be used on the output of \code{coalesc()} or \code{forward()} functions.
}
\usage{
plot.comm(x, type = "trait", seltrait = 1, main = NULL)
}
\arguments{
 \item{x}{
A list including the species pool composition (\code{x$pool}) and the local community composition (\code{x$com})
}
 \item{type}{
If \code{"trait"}, the function displays density plots of trait distributions. If \code{"abund"}, it displays the relationship between local and regional abundances
}
 \item{seltrait}{
Index of the trait to be plotted (if multiple traits used in simulation)
}
 \item{main}{
A title to add on the plot
}
}
\details{
This generic function should be improved in the future for additional options of graphical display
}
\value{
Return two stacked density plots
}
\author{
François Munoz
}

\examples{
# Simulation of a neutral community including 100 individuals
J <- 500; theta <- 50; m <- 0.1;
comm1 <- coalesc(J, m, theta) 
plot.comm(comm1)

# Stabilizing habitat filtering around t = 0.5
comm2 <- coalesc(J, m, theta, filt = function(x) 0.5 - abs(0.5 - x))
plot.comm(comm2)
}

\keyword{trait distribution}
\keyword{species abundances}
\keyword{regional pool}
\keyword{local community}
