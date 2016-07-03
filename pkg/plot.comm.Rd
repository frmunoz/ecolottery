\name{plot.comm}
\alias{plot.comm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting regional and local trait distributions
}
\description{
Density plots of the trait distributions in the regional pool and in a local community. To be used on the output of coalesc or forward functions 
}
\usage{
plot.comm(x, seltrait=1)
}
\arguments{
 \item{x}{
A list including the species pool composition ($pool) and the local community composition ($com)
}
 \item{seltrait}{
Index of the trait to be plotted (if multiple traits used in simulation).
}
}
\details{
This generic function should be improved in the future for additional options of graphical display.
}
\value{
}
\references{
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
\keyword{regional pool}
\keyword{local community}
