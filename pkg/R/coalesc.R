coalesc <- function(J, theta, m = 1, filt = NULL, pool = NULL, Jpool = 50*J) {
  #Create the regional pool if not provided
  if (is.null(pool)) {
    AA <- Jpool # Total number of individuals in the pool
    A2 <- 1:Jpool  # Labels of individuals
    S2 <- array(0, c(AA, 1))  # Species labels
    C2 <- array(0, c(AA, 1))  # Trait
    aa <- 1  # ?
    ss <- 0  # ?
    Y <- runif(AA)  # Generate a vector to determine species
    
    R2 <- theta / (theta + (1:AA) - 1)  # Vector to determine species
    mia <- which(Y <= R2)  # Get all individuals with different species
    # All individuals that are reassigned to previous secies
    mib <- which(Y > R2)
    
    S2[mia] <- 1:length(mia)  # Set species number
    C2[mia] <- runif(length(mia))  # Compute Trait
    
    # For all individuals without an assigned species
    for (j in 1:length(mib)) {
      # Select randomly a previously assigned individual
      jj <- max(round(runif(1) * mib[j] - 1), 1)
      A2[mib[j]] <- A2[jj]  # ??? Replace value in individual list
      S2[mib[j]] <- S2[jj]  # Assign species of previously assigned individual
      C2[mib[j]] <- C2[jj]  # Assign species trait
    }
    pool <- cbind(A2, S2, C2)
  } else if (ncol(pool) < 2) {
    stop("The regional pool is misdefined (at least two columns required)")
  } else if (ncol(pool) == 2) {
    warning("No trait information provided in the regional pool")
  }
  
  if (!is.null(filt)) {
    filt1 <- function(x) filt(x)
  } else {
    filt1 <- function(x) sapply(x, function(x) 1)
  }
  
  
  if (!is.null(filt) & ncol(pool) == 2) {
    pool[, 3] <- rep(1, nrow(pool))
  }
  
  A <- array(0, c(J, 1))
  S <- array(0, c(J, 1))
  C <- array(0, c(J, 1))
  
  if (m < 1) {
    I <- m *(J - 1) / (1 - m)
    X <- runif(J)
    R1 <- I / (I + (1:J) - 1)
    mi1 <- which(X <= R1)
    mi2 <- which(X > R1)
  } else {
    mi1 = 1:J; mi2 = c()
  }
  
  AA <- length(mi1)
  
  migrants <- sample(1:nrow(pool), AA, prob = filt1(pool[, 3]))
  A[mi1] <- pool[migrants, 1]
  S[mi1] <- pool[migrants, 2]
  C[mi1] <- pool[migrants, 3]
  
  for (j in 1:length(mi2)) {
    if (j > 1) {
      
      if (A[mi2[j]] != 0) {
        stop("Error in the assignation of ancestors")
      }
      
      jj <- sample(c(mi1[mi1 < mi2[j]], mi2[1:(j - 1)]), 1)
    }
    else {
      jj <- sample(mi1[mi1 < mi2[j]], 1)
    }
    
    A[mi2[j]] <- A[jj]
    S[mi2[j]] <- S[jj]
    C[mi2[j]] <- C[jj]
  }
  
  res <- cbind(A, S, C)
  
  return(res)

}
