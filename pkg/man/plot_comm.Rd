\name{plot_comm}
\alias{plot_comm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plotting trait and species distributions of simulated communities
}
\description{
  Graphical function to used on the output of \code{coalesc()} or
  \code{forward()} functions. It can show the links between regional and
  local trait/abundance distributions, or local species and rank
  abundance distributions.
}
\usage{
plot_comm(x, type = "trait", seltrait = 1, main = NULL)
}
\arguments{
 \item{x}{
 a list including the species pool composition (\code{x$pool}) and the local
 community composition (\code{x$com}). For example, x may be the output of
 \code{coalesc()} or \code{forward()} functions.
}
 \item{type}{
   \itemize{
    \item if \code{type = "trait"}, the function displays density plots of trait
          distributions.
    \item if \code{type = "locreg"},it displays the relationship between local and regional abundances.
    \item if \code{type = "sad"}, it displays the empirical percentiles and theoretical values of species
          abundance distribution (SAD) in the community.
    \item if \code{type = "rad"}, it displays the empirical percentiles and theoretical values of rank
          abundance distribution (RAD) in the community.
   }
}
 \item{seltrait}{
  index of the trait to be plotted following community data.frame (if multiple
  traits used in simulation).
}
 \item{main}{
  an overall title for the plot.
}
}
\details{
  If \code{type = "trait"}, the function provides density plots of the trait or
  abundance distributions in the regional pool and in a local community.
  If \code{type = "locreg"}, it displays the relationship between
  regional and local species relative abundances.
  If \code{type = "sad"} or \code{"rad"}, it shows the percentile plots provided by 
  functions in the \code{"sads"} package (see \code{\link[sads]{ppsad}}). The reference red line corresponds to log-series 
  for \code{"sad"} and to geometric series for \code{"rad"}.
  
  By default \code{type = "trait"}.
  
  To be used on the output of \code{coalesc()} or \code{forward()} functions.
}
\value{
Return two stacked \code{\link[=ggplot2-package]{ggplot2}} density plots if
\code{type = "trait"} and biplots otherwise.
}
\author{
F. Munoz; P. Denelle
}

\examples{
# Simulation of a neutral community including 100 individuals
J <- 500; theta <- 50; m <- 0.1;
comm1 <- coalesc(J, m, theta) 
plot_comm(comm1)
plot_comm(comm1, type = "locreg")
plot_comm(comm1, type = "sad")

# Stabilizing habitat filtering around t = 0.5
filt_gaussian <- function(x) exp(-(x - 0.5)^2/(2*0.1^2))
comm2 <- coalesc(J, m, theta, filt = filt_gaussian)
plot_comm(comm2)
plot_comm(comm2, type = "locreg")
plot_comm(comm1, type = "sad")
}

\keyword{trait distribution}
\keyword{species abundances}
\keyword{regional pool}
\keyword{local community}
