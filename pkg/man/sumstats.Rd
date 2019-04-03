\name{sumstats}
\alias{sumstats}
\title{
Predifined summary statistic function for estimating neutral and non-neutral parameters of community assembly using Approximate Bayesian Computation (ABC) implemented in coalesc_abc2. 
}
\description{
Adaptable statistic function to be provided into \code{coalesc_abc2()} to estimate parameters of neutral migration-drift dynamics and parameters of environmental filtering using ABC and adapted to the format imput of the local community. These return a set of diversity indexes which can be based on species local abundances, functional traits or both. 
}
\usage{
  sumstats(com, multi="single", traits = NULL, type = "mix", n = 4, m = 4)
  }
  
\arguments{
  \item{com}{
    the observed community composition. If \code{multi = "single"} (default),
    should be a matrix or data.frame of individuals on rows with their
    individual id (first column), species id (second column) and 
    trait information (optional third column)
      }
  \item{multi}{
    structure of the community input:
      \itemize{
        \item if \code{multi = "single"}, \code{comm.obs} contains a single
        community. 
        
        \item if \code{multi = "tab"}, the user provides a site-species matrix
        (\strong{sites in rows and species in columns})
        
        \item if \code{multi = "seqcom"}, \code{comm.obs} contains a list of
        communities (not yet)
      }
    }
  \item{traits}{
    trait values (one or several traits in column) for the species present in \code{com}.
  }
  \item{type}{
  Type of diversity indices to be included in the summary statistics:
      \itemize{
        \item if \code{type = "taxo"}, summary statistics are the \code{n} first Hill 
        diversity numbers 
        % given by: N_a = \sum{i=1}^S (p_i^a)^{1/(1-a)} with 0 < a < n 
        
        \item if \code{type = "func"}, summary statistics are the 4 first moments of 
        the trait distribution, corresponding to the community weighted mean, variance, 
        skewness and kurtosis.  
        
        \item if \code{type = "mix"} (default),summary staitstics include the first 4 moments
        of the trait distribution (cwm, cwv, cws, cwk) as well as the first 4 Hill diversity 
        numbers. 
      }
  }
  \item{n}{
    number of Hill's diversity numbers to compute when \code{type} is taxo or mix. 
    Default is \code{n=4}. 
  }
  
  \item{m}{
    number of trait distribution moments to compute when \code{type} is func or mix. 
    Default is \code{m=4}. 
  }
}

\details{
  \code{sumstats} provides summary statistics based on user preference and the format of community
  data to be used in \code{coalesc_abc2} for investigating neutral or niche-based dynamics in 
  community assembly using ABC. 
  
}
\value{
A vector of summary statistics to be used to compare observed and simulated community composition in the ABC estimation performed by \code{coalesc_abc2}. 

For community information of the type \code{multi = "single"} this vector is comprised of n Hill numbers (when \code{type} is "taxo", the m first moments of the trait distribution when \code{type} is "func" or the m first moments of the trait distribution followed by the n Hill numbers by  when \code{type} is "mix" of the single community. 

For community information as a species by site matrix or data.frame with X sites, this vector is comprised of n Hill numbers of the first community, followed by the n Hill numbers for the second and so on until X (when \code{type} is "taxo"). Likewise (when \code{type} is "func", with m first moments of the trait distributions of the X communities. When \code{type} is "mix", the vector gives the m first moments of the trait distribution for every community followed by their n Hill numbers. 
}

\references{
  Ref? 
  
}
\author{
  E.Barthelemy
}
\examples{
  \dontrun{
    examples? 
   } 
}

\keyword{summary statistics}
\keyword{Approximate Bayesian Computation}
\keyword{niche-based dynamics}
\keyword{neutral dynamics}
