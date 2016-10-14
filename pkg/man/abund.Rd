\name{abund}
\alias{abund}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of local and regional abundances and relative abundances
}
\description{
Calculate the abundances and relative abundances of species in simulated communities and in the corresponding species pools
}
\usage{
abund(x)
}
\arguments{
\item{x}{
A list including the species pool composition ($pool) and the local community composition ($com)
}
}
\details{
This generic function should be improved in the future for additional options.
}
\value{
\item{reg}{
Species abundances and relative abundances in the regional pool
}
\item{loc}{
Species abundances and relative abundances in the local community
}
}
\author{
François Munoz
}

\examples{
J <- 500; theta <- 50; m <- 0.05;
comm1a <- coalesc(J, m, theta) 
abund1a <- abund(comm1a) 
# Log-series distribution of regional abundances
require(vegan)
plot(fisherfit(abund1a$reg$ab))

# Relationship between local and regional abundances
par(mfrow=c(1,2))
plot(abund1a$reg$abrel[abund1a$loc$sp],abund1a$loc$abrel,main="m = 0.1",xlab="Regional abundance", ylab="Local abundance",log="xy"); abline(0,1)
# With higher immigration rate
m <- 0.95; comm1b <- coalesc(J, m, theta) 
abund1b <- abund(comm1b) 
plot(abund1b$reg$abrel[abund1b$loc$sp],abund1b$loc$abrel,main="m = 0.1",xlab="Regional abundance", ylab="Local abundance",log="xy"); abline(0,1)
}

\keyword{species abundance distribution}
\keyword{regional pool}
\keyword{local community}
