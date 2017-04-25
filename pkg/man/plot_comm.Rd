\name{plot_comm}
\alias{plot_comm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Regional vs. Local trait distributions of abundances
}
\description{
  Graphical function to used on the output of \code{coalesc()} or
  \code{forward()} functions.It aims at plotting links between regional and
  local trait/abundance distributions.
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
    \item if \code{type = "abund"}, it displays the relationship between local
          and regional abundances.
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
  If \code{type = "abund"}, the function displays the relationship between
  regional and local species relative abundances.
  By default \code{type = "trait"}.
  To be used on the output of \code{coalesc()} or \code{forward()} functions.
}
\value{
Return two stacked \code{\link[=ggplot2-package]{ggplot2}} density plots if
\code{type = "trait"} and a biplot if \code{type = "abund"}.
}
\author{
F. Munoz; P. Denelle
}

\examples{
# Simulation of a neutral community including 100 individuals
J <- 500; theta <- 50; m <- 0.1;
comm1 <- coalesc(J, m, theta) 
plot_comm(comm1)
plot_comm(comm1, type = "abund")

# Stabilizing habitat filtering around t = 0.5
comm2 <- coalesc(J, m, theta, filt = function(x) 0.5 - abs(0.5 - x))
plot_comm(comm2)
plot_comm(comm2, type = "abund")
}

\keyword{trait distribution}
\keyword{species abundances}
\keyword{regional pool}
\keyword{local community}
