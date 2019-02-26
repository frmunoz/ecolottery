coalesc <- function(J, m = 1, theta = NULL, filt = NULL, add = F,  var.add =NULL, pool = NULL,
                    traits = NULL, Jpool = 50*J, verbose = FALSE) {
  
  # Check parameters
  if (is.null(theta) & is.null(pool)) {
    stop("You must provide either regional pool composition or a theta value")
  }
  
  if((add & is.null(var.add)) | (!add & !is.null(var.add))) {
    warning("No additional variables are passed to filt")
  }
    
  if(length(m)>1 | length(theta)>1) {
    stop("m and theta cannot be vectors of length greater than 1")
  }
  
  if (m < 0 | m > 1) {
    stop("The migration parameter takes values between 0 and 1")
  }
  
  if (!is.null(theta)) {
    if (theta <= 0) {
      stop("The theta parameter must be positive")
    } else if (theta > 0 & !is.null(pool)) {
      if (verbose) {
        warning("Both a theta value and a regional pool provided, discarding ",
                "theta")
      }
    }
  }
  
  if (J <= 0) {
    stop("J must be positive")
  }
  
  if (is.null(traits) & (is.null(pool) | NCOL(pool) < 3)) {
    if (verbose) warning("No trait information provided in the regional pool")
  }
  
  if (!is.null(traits) & is.null(colnames(traits))) {
    colnames(traits) <- paste("tra", 1:ncol(traits), sep = "")
  }
  if (!is.null(pool) & is.null(colnames(pool))) {
    if (ncol(pool) > 2) {
      colnames(pool) <- c("ind", "sp", paste("tra", 1:(ncol(pool) - 2),
                                             sep = ""))
    }
  }
  
  ## Create the regional pool if not provided
  
  if (is.null(pool)) {
    # If m = 1 and no environmental filtering, directly simulates a sample from
    # the pool of size J
    if (m == 1 & is.null(filt)) pool_size <- J else  pool_size <- Jpool
    ind_pool_lab <- 1:pool_size  # Labels of individuals
    ind_pool_sp <- array(NA, c(pool_size, 1))  # Species labels
    
    if (is.null(traits)) {
      ind_pool_traits <- array(NA, c(pool_size, 1))
    } else {
      ind_pool_traits <- array(NA, c(pool_size, ncol(traits)))
    }
    
    ## Determining species identity in the pool
    Y <- runif(pool_size)  # Generate a vector to determine species
      
    # Probability that new species are sampled in the regional pool
    R_pool <- theta / (theta + (1:pool_size) - 1) 
    # Sampled individuals belonging to new species
    assign_pool <- which(Y <= R_pool)
    # Other individuals are assigned to already sampled species
    unassign_pool <- which(Y > R_pool)
    
    ind_pool_sp[assign_pool] <- 1:length(assign_pool)  # Set species number
    
    if (is.null(traits)) {
      # Assign random trait values
      ind_pool_traits[assign_pool,1] <- runif(length(assign_pool))  
    } else {
      ind_pool_traits[assign_pool,] <- traits[1:length(assign_pool),]
    }

    # For all individuals without an assigned species
    if (length(unassign_pool) > 0) for (j in 1:length(unassign_pool)) {
      # Select randomly a previously assigned individual
      existing_sp <- sample.int(unassign_pool[j] - 1, 1)
      # Assign species of previously assigned individual
      ind_pool_sp[unassign_pool[j]] <- ind_pool_sp[existing_sp] 
      # Assign species trait
      ind_pool_traits[unassign_pool[j],] <- ind_pool_traits[existing_sp,]
    }
    
    if (!is.null(traits)) {
      colnames(ind_pool_traits) <- colnames(traits)
    } else {
      colnames(ind_pool_traits) <- paste("tra", 1:ncol(ind_pool_traits),
                                         sep = "")
    }
    
    if (m == 1 & is.null(filt)) {
      return(list(pool = data.frame(ind = ind_pool_lab, sp = ind_pool_sp,
                                    ind_pool_traits), call = match.call()))
    }
  } else { 
    
    if (ncol(pool) < 2) {
      stop("The regional pool is misdefined (two columns required)")
    }
    
    # Handling of factors
    if (is.factor(pool[,1])) pool[,1] <- as.character(pool[,1])
    if (is.factor(pool[,2])) pool[,2] <- as.character(pool[,2])
      
    ind_pool_lab <- pool[,1]; ind_pool_sp <- pool[,2]
    
    if (ncol(pool) >= 3)  {
      ind_pool_traits <- data.frame(pool[,-(1:2)])
      if(any(is.na(ind_pool_traits)))
        stop("There should not be any NA in trait values of the individuals of the pool")
      if(!is.null(traits))  
        warning("Trait information already in 'pool', the 'traits' input is ignored")
      colnames(ind_pool_traits) <- colnames(pool)[-(1:2)]
    } else {
      # Generation of trait values if not provided by the user
      if (is.null(traits)) {
        traits <- data.frame("tra" = runif(nlevels(as.factor(pool[,2]))))
        rownames(traits) <- levels(as.factor(pool[,2]))
      } else {
        if(!is.null(rownames(traits)) & any(!ind_pool_sp %in% rownames(traits))) 
          stop("Species names in traits must match those in pool")
      }
      
      ind_pool_traits <- data.frame(traits[ind_pool_sp,1:ncol(traits)], stringsAsFactors = F)
      
      if (!is.null(traits)) {
        colnames(ind_pool_traits) <- colnames(traits)
      } else {
        colnames(ind_pool_traits) <- paste("tra", 1:ncol(ind_pool_traits),
                                           sep = "")
      }
    }
  }
  
  pool <- data.frame(ind = ind_pool_lab, sp = ind_pool_sp, ind_pool_traits)
  
  for (i in which(unlist(lapply(pool,is.factor)))) {
    pool[,i] <- as.character(pool[,i])
  }
  
  ## Define environmental filter
        
  if (!is.null(filt)) {
    if(!add)
    { #env_filter <- function(x) t(apply(x, 1, function(y) filt(y)))
      env_filter <- function(x) unlist(sapply(1:nrow(x), function(i) filt(x[i,])))
      
    } else 
    {
      #env_filter <- function(x, var.add) t(apply(x, 1, function(y) filt(y, var.add)))
      env_filter <- function(x) unlist(sapply(1:nrow(x), function(i) filt(x[i,], var.add)))
      
    }
  } else {
     #env_filter <- function(x) t(apply(x, 1, function(y) 1))
     env_filter <- function(x) unlist(sapply(1:nrow(x), function(i) 1))
    
  }
  
  # Environmental filter should not provide negative values
  
  if (!add) 
  {
    prob <- env_filter(ind_pool_traits)
  } else prob <- env_filter(ind_pool_traits, var.add)
  
  if (any(prob < 0)) {
      if (verbose) {
        warning("Negative weights yielded by filtering function are set to ",
                "0.\n", "Maybe better defining the filtering function in a ",
                "different way!")
      }
      prob[prob < 0] <- 0
  }  
  
  if(all(prob == 0)) {
    if (verbose) {
      warning("Your filtering function does not allow any immigrant to enter ",
              "the community")
    }
    return(list(com = array(NA, c(J,3)), pool = pool, call = match.call()))
  }

  ## Community generation
  
  ind_com_lab <- array(NA, J)
  ind_com_sp <- array(NA, J)
  ind_com_traits <- matrix(NA, nrow = J, ncol = ncol(ind_pool_traits))      
  
  # If migration rate = 0, ecological fixation
  # If J = 1, the community is made of an individual drawn from the pool
  if (m == 0 | J == 1) {
    
    migrant <- sample(1:nrow(pool), 1, prob = prob)
    
    ind_com_traits <- t(apply(ind_com_traits, 1,
                              function(x) ind_pool_traits[migrant, ]))
    
    com <- data.frame(ind = rep(pool[migrant, 1], J),
                      sp = rep(pool[migrant, 2], J),
                      ind_com_traits)
    return(list(com = com, pool = pool, call = match.call()))
  }
    
  # If migration rate > 0 and more than 1 individual in community
  if (m < 1 & J > 1) {
    I <- m * (J - 1) / (1 - m) # Number of available immigrants
    X <- runif(J)
    # Probability that new migrant ancestor arrives in the community
    R_com <- I / (I + (1:J) - 1)
    assign_com <- which(X <= R_com)
    unassign_com <- which(X > R_com)
  } else {
    # migration rate set at its maximum, unlimited regional pool.
    # All new individuals taken from the regional pool
    assign_com <- 1:J
    unassign_com <- NULL
  }
  
  # Number of individuals taken from the regional pool
  com_species <- length(assign_com) 
 
  if(com_species > sum(prob > 0)) {
    warning("Sampling with replacement from the pool")
    migrants <- sample(1:nrow(pool), com_species, replace = T, prob = prob)
  } else 
    migrants <- sample(1:nrow(pool), com_species, prob = prob)
  
  # Assign individuals
  ind_com_lab[assign_com] <- pool[migrants[1:com_species], 1] 
  
  ind_com_sp[assign_com] <- pool[migrants[1:com_species], 2] # Assign species
 
  # Assign traits
  ind_com_traits[assign_com,] <- apply(ind_pool_traits, 2,
                                       function(y) y[migrants[1:com_species]]) 
  colnames(ind_com_traits) <- colnames(ind_pool_traits)
 
  if (!is.null(unassign_com) & length(unassign_com) > 0) {
    for (j in 1:length(unassign_com)) {
    existing_sp <- sample.int(unassign_com[j] - 1, 1)
    ind_com_lab[unassign_com[j]] <- ind_com_lab[existing_sp]
    ind_com_sp[unassign_com[j]] <- ind_com_sp[existing_sp]
    ind_com_traits[unassign_com[j],] <- ind_com_traits[existing_sp,]
    }
  }
  
  # For debug only
  if (any(is.na(ind_com_lab))) {
    stop("NA in simulated community")
  }
    
  com <- data.frame(ind = ind_com_lab, sp = ind_com_sp, ind_com_traits, stringsAsFactors = F)
  
  if (m == 1 & is.null(filt)) {
    return(list(pool = com, call = match.call()))
  } else {
    return(list(com = as.data.frame(com), pool = as.data.frame(pool), call = match.call()))
  }
}
