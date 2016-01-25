coalesc <- function(J, theta, m = 1, filt = NULL, pool = NULL, Jpool = 50*J) {
  #Create the regional pool if not provided
  if (is.null(pool)) {
    
    pool_size <- Jpool # Total number of individuals in the pool
    ind_pool_lab <- 1:Jpool  # Labels of individuals
    sp_pool_lab <- array(0, c(pool_size, 1))  # Species labels
    sp_trait <- array(0, c(pool_size, 1))  # Trait
    Y <- runif(pool_size)  # Generate a vector to determine species
    
    # Vector to determine species
    R_pool <- theta / (theta + (1:pool_size) - 1)
    assign_pool <- which(Y <= R_pool)  # Get all individuals with different species
    
    # All individuals that are reassigned to previous secies
    unassign_pool <- which(Y > R_pool)
    
    sp_pool_lab[assign_pool] <- 1:length(assign_pool)  # Set species number
    sp_trait[assign_pool] <- runif(length(assign_pool))  # Compute Trait

    # For all individuals without an assigned species
    for (j in 1:length(unassign_pool)) {
      # Select randomly a previously assigned individual
      jj <- max(round(runif(1) * unassign_pool[j] - 1), 1)
      
      ind_pool_lab[unassign_pool[j]] <- ind_pool_lab[jj]  # ??? Replace value in individual list
      sp_pool_lab[unassign_pool[j]] <- sp_pool_lab[jj]  # Assign species of previously assigned individual
      sp_trait[unassign_pool[j]] <- sp_trait[jj]  # Assign species trait

    }
    pool <- cbind(ind_pool_lab, sp_pool_lab, sp_trait)
    
  } else if (ncol(pool) < 2) {
    stop("The regional pool is misdefined (at least two columns required)")
  } else if (ncol(pool) == 2) {
    warning("No trait information provided in the regional pool")
  }
  
  
  # Define environmental filter
  if (!is.null(filt)) {
    env_filter <- function(x) filt(x)
  } else {
    env_filter <- function(x) sapply(x, function(x) 1)
  }
  
  # No environmental filter and no trait defined, all traits == 1

  if (!is.null(filt) & ncol(pool) == 2) {
    pool[, 3] <- rep(1, nrow(pool)) 
  }
  
  # Community Array
  ind_com_lab <- array(0, c(J, 1))
  sp_com_lab <- array(0, c(J, 1))
  sp_com_trait <- array(0, c(J, 1))

  # If migration rate < 1
  if (m < 1) {
    I <- m * (J - 1) / (1 - m) # Number of available immigrants
    X <- runif(J)
    # Probability that new migrant ancestor arrives in the community
    R_com <- I / (I + (1:J) - 1)
    assign_com <- which(X <= R_com)
    unassign_com <- which(X > R_com)
  } else {
    # migration rate set at its maximum, unlimited regional pool.
    # All new individuals taken from the regional pool
    assign_com = 1:J
    unassign_com = NULL
  }
  
  com_species <- length(assign_com) # Number of individuals taken from the regional pool
  
  migrants <- sample(1:nrow(pool), com_species, prob = env_filter(pool[, 3]))
  
  ind_com_lab[assign_com] <- pool[migrants, 1] # Assign individuals
  
  sp_com_lab[assign_com] <- pool[migrants, 2] # Assign species
  
  sp_com_trait[assign_com] <- pool[migrants, 3] # Assign traits
  
  if (!is.null(unassign_com)) {
    for (j in 1:length(unassign_com)) {
      if (j > 1) {
        
        if (ind_com_lab[unassign_com[j]] != 0) {
          stop("Error in the assignation of ancestors")
        }
        
        jj <- sample(c(assign_com[assign_com < unassign_com[j]], unassign_com[1:(j - 1)]), 1)
      }
      else {
        jj <- sample(assign_com[assign_com < unassign_com[j]], 1)
      }
      
      ind_com_lab[unassign_com[j]] <- ind_com_lab[jj]
      sp_com_lab[unassign_com[j]] <- sp_com_lab[jj]
      sp_com_trait[unassign_com[j]] <- sp_com_trait[jj]
    }
  }
    
  com <- cbind(ind_com_lab, sp_com_lab, sp_com_trait)
  
  return(com)
}

