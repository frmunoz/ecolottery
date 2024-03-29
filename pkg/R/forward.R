# Function to compute forward simulation of community dynamics with (eventually)
# environmental filtering
forward <- function(initial, m = 1, theta = NULL, d = 1, gens = 150, keep = FALSE,
                    pool = NULL, traits = NULL, filt = NULL, filt.vect = F,
                    limit.sim = NULL, limit.intra = F, par.limit = 0.1, coeff.lim.sim = 1, 
                    type.filt = "immig", type.limit = "death", m.replace = T, add = F, var.add = NULL,
                    prob.death = NULL, method.dist = "euclidean", checks = T, plot_gens = FALSE) {
  # The function will stop if niche - based dynamics is requested, but trait
  # information is missing in the local community
  # For strictly neutral communities, a vector of species names is enough for
  # the initial community
  
  if(checks)
  {
    if (is.null(theta) & is.null(pool)) {
      stop("You must provide either regional pool composition or a theta value",
           call. = FALSE)
    }
    
    # Need to define a function for environmental filtering
    if(!is.null(filt)) if(!is.function(filt)) filt <- NULL
    if(is.null(filt)) filt.vect <- F
    
    if((add & is.null(var.add)) | (!add & !is.null(var.add))) {
      warning("No additional variables are passed to filt", call. = FALSE)
    }
    
    if(length(m)>1 | length(theta)>1) {
      stop("m and theta cannot be vectors of length greater than 1",
           call. = FALSE)
    }
    
    if (m < 0 | m > 1) {
      stop("The migration parameter takes values between 0 and 1", call. = FALSE)
    }
    
    if (!is.null(theta)) {
      if (theta <= 0) {
        stop("The theta parameter must be positive", call. = FALSE)
      } else if (theta > 0 & !is.null(pool)) {
          warning("Both a theta value and a regional pool provided, discarding ",
                  "theta", call. = FALSE)
      }
    }
   
    # When composition is given as a vector of species names
    if(is.factor(initial))
      initial <- as.character(initial)
    if(is.factor(pool))
      pool <- as.character(pool)
    
    # Standard formatting: first column is individual id, second is species id
    if(is.vector(initial))
      initial <- data.frame(ind = 1:length(initial), sp = initial, stringsAsFactors = FALSE)
    if(is.vector(pool)) {
      pool <- data.frame(ind = paste("pool",1:length(pool),sep="."), sp = pool, stringsAsFactors = FALSE)
      rownames(pool) <- pool$ind
    }
    
    if (!is.null(pool) & is.null(colnames(pool))) {
      if (ncol(pool) > 2) {
        colnames(pool) <- c("ind", "sp", paste("tra", 1:(ncol(pool) - 2),
                                               sep = ""))
      }
    }
    
    if (is.null(traits) & (is.null(pool) | NCOL(pool) < 3)) {
      warning("No trait information provided in the regional pool",
                           call. = FALSE)
    }
    
    if (!is.null(traits) & is.null(colnames(traits))) {
      colnames(traits) <- paste("tra", 1:ncol(traits), sep = "")
    }
    
    if (!is.null(pool)) if (ncol(pool) < 2) {
      stop("The regional pool is misdefined (at least two columns ",
           "required when a matrix or data frame is provided)", call. = FALSE)
    }
    
    #if(any(initial$ind %in% pool$ind))
    #  warning("It is advisable not to have the same individual ids in pool and initial")
    
    # For back-compatibility
    if(is.logical(limit.sim)) 
    {
      if(!limit.sim) limit.sim <- NULL else limit.sim <- gauss_limit
    }
    
    if(!is.null(limit.sim) & !is.function(limit.sim))
      stop("limit.sim must either be NULL or a function defining limiting similarity")
    
    if (!is.numeric(d) | d < 0) {
      stop("Number of individuals that die in each time step must be a positive",
           " number.", call. = FALSE)
    }
    
    if (!is.numeric(gens) | gens <= 0) {
      stop("Number of generations must be a positive number.", call. = FALSE)
    }
    
    if (!is.logical(keep)) {
      stop("keep parameter must be a boolean.")
    }
    
    if (!is.numeric(coeff.lim.sim)) {
      stop("coeff.lim.sim parameter must be numeric.", call. = FALSE)
    }
    
    if (!is.null(filt) & (is.null(type.filt) | !any(type.filt%in%c("immig","death","loc.recr"))))
      stop("Type of environmental filtering should be immig, death and/or ",
           "loc.recr", call. = FALSE)
    
    if (!is.null(limit.sim) & (is.null(type.limit) | !any(type.limit%in%c("immig","death","loc.recr"))))
      stop("Type of limiting similarity should be immig, death and/or loc.recr",
           call. = FALSE)
    
    if (!is.null(filt) & !is.function(filt)) {
      stop("filt() must be a function.", call. = FALSE)
    }    
    
    if(add & is.null(var.add)) 
      warning("No additional variables are passed to filt", call. = FALSE)
    
    if ((method.dist %in% c("euclidean", "maximum", "manhattan", "canberra",
                            "binary", "minkowski")) == FALSE) {
      stop("Provided distance does not exist. See stats::dist function for help.",
           call. = FALSE)
    }
    
    if (!is.logical(plot_gens)) {
      stop("plot_gens parameter must be a boolean.", call. = FALSE)
    }
    
    # If no limiting similarity nor environmental filter: neutral community dynamics 
    if (is.null(limit.sim) & is.null(filt)) {
      message("Simulation of a neutral community")
    }
    
    if (!is.null(traits) & is.null(colnames(traits))) {
      colnames(traits) <- paste("tra", 1:ncol(traits), sep = "")
    }
    
    # Check and complete trait values for trait-based dynamics
    if (!is.null(limit.sim) | !is.null(filt)) {
      # Trait information in pool
      if(!is.null(pool)) if(ncol(pool) == 2 &  !is.null(traits)) {
        # Assign trait values to the individuals of the pool
        pool[,3:(2+ncol(traits))] <- traits[pool[,2],]
        colnames(pool)[-(1:2)] <- colnames(traits)
        if(any(is.na(pool[,-(1:2)]))) 
          stop("Mismatch of species names between pool and traits", call. = FALSE)
      }
      # Trait information in initial
      if(ncol(initial) == 2 &  !is.null(traits)) {
        # Assign trait values to the individuals of the initial community
        initial[,3:(2+ncol(traits))] <- traits[initial[,2],]
        colnames(initial)[-(1:2)] <- colnames(traits)
        if(any(is.na(initial[,-(1:2)]))) 
          stop("Mismatch of species names between initial and traits", call. = FALSE)
      }
      # Generation of random trait values
      if(is.null(traits) & ncol(initial) < 3)
        stop("Trait values are required in initial community")
      if(!is.null(pool)) if(is.null(traits) & ncol(pool) < 3) {
        for(i in 1:(ncol(initial)-2))
          pool[, i+2] <- runif(nrow(pool))
        colnames(pool)[-(1:2)] <- colnames(initial)[-(1:2)]
        message("Random (uniform) trait values attributed to individuals of ",
                "the regional pool")
      }
    }
    
    if (!is.null(pool) & is.null(colnames(pool))) {
      if (ncol(pool) > 2) {
        colnames(pool) <- c("ind", "sp", paste("tra", 1:(ncol(pool) - 2),
                                               sep = ""))
      }
    }
    
    J <- nrow(initial)
   
    if(is.null(filt) & !is.null(pool)) if(ncol(pool)==2 & ncol(initial)>2)
    {
      # TO BE IMPROVED
      # No need to keep trait information for neutral dynamics
      initial <- initial[,1:2]
    }
    
    if (J < d) stop("The number of dead individuals per time step ",
                    "cannot be greater than community size", call. = FALSE)
    
    if ((!is.null(limit.sim) | !is.null(filt))) if(any(is.na(initial[, 3]))) {
      stop("Trait information must be provided in initial community ",
           "composition for niche-based dynamics", call. = FALSE)
    }
  }
  
  new.index <- 0
  
  ## Forward simulation with community
  
  # Begins with the initial community
  next_comm <- initial
  
  # Richness of initial community is not included
  sp_t <- c()
  
  ind_t <- c()	
  dist.t <- c()
  
  if (keep) {  # If the user asked to keep all the communities at each timestep
    comm_through_time <- c()
  }
  
  # Simulate the community for the given number of generations
  for (i in 1:gens) {
    if (keep) {
      # Store the community at time i
      comm_through_time[[i]] <- next_comm
    }  
    
    # Simulate community dynamics
    next_comm <- pick(next_comm, d = d, m = m, theta = theta, pool = pool,
                      prob.death = prob.death, filt = filt, filt.vect = filt.vect,
                      limit.sim = limit.sim, limit.intra = limit.intra,
                      par.limit = par.limit, 
                      coeff.lim.sim = coeff.lim.sim, 
                      type.filt = type.filt, type.limit = type.limit, 
                      m.replace = m.replace,
                      add = add, var.add = var.add,
                      new.index = new.index,
                      method.dist = "euclidean")
    
    sp_t <- c(sp_t, length(unique(next_comm$com$sp)))
    ind_t <- c(ind_t, length(unique(next_comm$com$ind)))
    
    # Store average trait distance among coexisting individuals
    if (!is.null(limit.sim)) {
      dist.t <- c(dist.t, next_comm$dist.t)
    }
    new.index <- next_comm$new.index
    next_comm <- next_comm$com
  }
  
  if (plot_gens) { # Plotting number of individuals and species over generations
    
    uniq_df <- data.frame(gens = 1:gens, ind_t = ind_t, sp_t = sp_t,
                          stringsAsFactors = FALSE)
    
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      
      # Plot the number of individuals through all the generations
      plot_individuals <- ggplot(uniq_df, aes_string("gens", "ind_t")) +
        geom_line() +
        geom_line(aes_string("gens", "ind_t"),
                  size = 1) +
        labs(x = "Number of generations",
             y = "Number of distinct ancestors")
      
      # Plot the number of species through all the generations
      plot_species <- ggplot(uniq_df,
                             aes_string("gens", "sp_t")) +
        geom_line(size = 1) +
        labs(x = "Number of generations",
             y = "Number of species")
      
      print(plot_individuals)
      print(plot_species)
    }
  }
  
  if (!is.null(limit.sim))
  {
    if (keep) return(list(com_t = comm_through_time,
                          com = comm_through_time[[length(comm_through_time)]],
                          sp_t = sp_t,
                          dist.t = dist.t,
                          pool = pool,
                          call = match.call()))
    else return(list(com = next_comm, sp_t = sp_t, 
                     dist.t = dist.t, pool = pool, call = match.call()))
  } else
  {
    if (keep) return(list(com_t = comm_through_time,
                          sp_t = sp_t,
                          pool = pool,
                          call = match.call()))
    else return(list(com = next_comm, sp_t = sp_t, pool = pool, call = match.call()))
  }
}

# Precise function to simulate a single timestep by picking an individual in
# the pool or make an individual mutate
pick <- function(com, d = 1, m = 1, theta = NULL, pool = NULL, prob.death = NULL,
                 filt = NULL, filt.vect = F, limit.sim = NULL, limit.intra = F, 
                 par.limit = 0.1, coeff.lim.sim = 1, 
                 type.filt = "immig", type.limit = "death", 
                 m.replace = T, add = F, var.add = NULL, new.index = 0,
                 method.dist = "euclidean") {
  
  
  J <- nrow(com)
  
  if (is.null(pool)) {
    # If no species pool specified, mutate an individual
    return(pick.mutate(com,
                       d = d,
                       mu = theta/(2*J),
                       new.index = new.index)) 
    
  } else {
    
    if((!is.null(filt) | !is.null(limit.sim)) & m > 0 & any(is.na(pool[,-(1:2)]))) {
      stop("With environmental filtering, NA trait values not allowed in ",
           "regional pool", call. = FALSE)
    }
    
    # If there is a species pool make an individual immigrates
    return(pick.immigrate(com, d = d, m = m, pool = pool,
                          prob.death = prob.death, filt = filt, filt.vect = filt.vect, 
                          limit.sim = limit.sim, limit.intra = limit.intra, 
                          par.limit = par.limit, coeff.lim.sim = coeff.lim.sim, 
                          type.filt = type.filt, type.limit = type.limit,
                          m.replace = m.replace, add = add, var.add = var.add,
                          method.dist = "euclidean"))
  }
}

# Return community with mutated individual (= new species)
pick.mutate <- function(com, d = 1, mu = 0, new.index = 0) {
  
  if (is.vector(com)) {
    # If community only defined by species names
    J <- length(com)
    
    com <- data.frame(ind = paste("ind", 1:J, sep = ""),
                      sp = as.character(com),
                      trait = rep(NA, J),
                      stringsAsFactors = FALSE)
    
  } else if (is.matrix(com) | is.data.frame(com)) {
    # If the community has defined traits
    J <- nrow(com)
    
    if (is.matrix(com)) {
      com <- as.data.frame(com, stringsAsFactors = FALSE)
      com[, 1] <- as.character(com[, 1])
      com[, 2] <- as.character(com[, 2])	  
    }
    
  } else {
    stop("pick.mutate: misdefined community composition", call. = FALSE)
  }
  
  ## Simulate the dynamics of the community
  
  # Number of individuals who die at this timestep
  died <- sample(J, d, replace = FALSE)
  
  # How many of the dead individuals are replaced by mutated individuals
  mutated <- runif(length(died)) < mu
  
  # Number of mutated individuals
  n_mutated <- sum(mutated)
  
  # Dead individuals who did not mutate
  dead_non_mutated <- sum(!mutated)
  
  # Replace dead non mutated individuals by individuals from other species
  com[died[!mutated], ] <- com[sample(1:nrow(com), dead_non_mutated,
                                      replace = TRUE), ]
  
  # When mutation occurs
  if (n_mutated > 0) {
    
    # Attribute new species name to individuals who mutate
    com[died[mutated], 2] <- paste("new.sp", new.index + (1:n_mutated),
                                     sep = "")
    
    #com[died[mutated], 3] <- rep(NA, n_mutated)  # No trait values
    # Default = trait value drawn from uniform distribution between 0 and 1
    com[died[mutated], 3] <- runif(n_mutated)
    
    # Number of new species which appeared (next one will be new.index + 1)
    new.index <- new.index + n_mutated
  }
  
  return(list(com = com, new.index = new.index))
}

# Function to return individuals who immigrated from the species pool
# limit.sim = limiting similarity; filt = habitat filtering function
pick.immigrate <- function(com, d = 1, m = 1, pool,
                           prob.death = NULL, filt= NULL, filt.vect = F,
                           limit.sim = NULL, limit.intra = F, par.limit = 0.1,
                           coeff.lim.sim = 1, type.filt = "immig", 
                           type.limit = "death", m.replace = T, add = F,
                           var.add = NULL, method.dist = "euclidean") {
  
  if (is.vector(com)) {
    # If community only defined by species names
    J <- length(com)
    
    com <- data.frame(ind = paste("ind", 1:J, sep = ""),
                      sp = as.character(com),
                      trait = rep(NA, J),
                      stringsAsFactors = FALSE)
    
  } else if (is.matrix(com) | is.data.frame(com)) {
    # If the community has defined traits
    J <- nrow(com)
    
    if (is.matrix(com)) {
      com <- as.data.frame(com, stringsAsFactors = FALSE)
      colnames(com)[1:2] <- c("ind", "sp")
    }
    com[, 1] <- as.character(com[, 1])
    com[, 2] <- as.character(com[, 2])
    
  } else {
    stop("pick.immigrate: misdefined community composition", call. = FALSE)
  }
  
  # Function defining habitat filtering according to trait value
  if(!filt.vect) {
    env_filter <- ifelse(!is.null(filt), ifelse(!add, 
                                                function(x) apply(x, 1, filt), 
                                                function(x, var.add) apply(x, 1, function(i) filt(i, var.add))), 
                         function(x) rep(1, nrow(x))) 
  } else {
    # If the filtering function can be directly applied to a vector of trait values
    env_filter <- filt
  }
  
  # Traits distances used to simulate limiting similarity
  if (!is.null(limit.sim)) {
    if(!("immig" %in% type.limit))
    {
      tr.dist <- as.matrix(dist(com[, -(1:2)], method = method.dist))
      sp.lab <- com[, 2]
      # Column and row names are individual labels
      colnames(tr.dist) <- com[, 1]
      rownames(tr.dist) <- com[, 1]
    } else {
      tr.dist <- as.matrix(dist(c(com[, -(1:2)], pool[, -(1:2)]), method = method.dist))
      sp.lab <- c(com[, 2], pool[, 2])
      colnames(tr.dist) <- c(com[, 1], pool[, 1])
      rownames(tr.dist) <- c(com[, 1], pool[, 1])
    }
    diag(tr.dist) <- NA
    
    # Should intraspecific competition be included
    if(!limit.intra)
    {
      # In this case no limiting similarity within species
      for(sp in unique(sp.lab))
        tr.dist[sp.lab==sp, sp.lab==sp] <- NA
    }
    # dist.t will display the average trait distance among species
    # for the whole community at each generation
    if (min(dim(tr.dist))>1) {
      dist.t <- mean(tr.dist[com[, 1], com[, 1]], na.rm = TRUE)
    } else dist.t <- NA
    
  } else {
    tr.dist <- NULL
  }
  
  # Initial community
  com.init <- com
  
  if (is.null(limit.sim) & is.null(filt)) {
    
    died <- sample(J, d, replace = FALSE)
    
    com <- com[-died, ]
    
    if (any(is.na(com[, 1]))) {
      stop("NA values in community composition (1)", call. = FALSE)
    }
    
  } else {
    
    if(!is.null(limit.sim))
    {
      # lim_sim_function indicates the influence of limiting similarity
      # coeff.lim.sim modulates the strength of limiting similarity
      lim_sim_function <- function(d) coeff.lim.sim * limit.sim(d)
    }
    
    # Vector of the individual probability of dying
    if (is.null(prob.death)) {
      prob.death <- rep(1, nrow(com))
    }
    
    # Influence of limiting similarity on mortality
    if("death" %in% type.limit & !is.null(limit.sim) & !is.null(tr.dist))
    {
      # Under limiting similarity, mortality increases when an individual is more
      # similar to other resident individuals
      # For each species: compute death probability depending on limiting
      # similarity plus a baseline individual death probability
      
      if(!is.null(limit.sim)) if(length(formals(limit.sim))>1)
      {
        prob.death <- prob.death + apply(tr.dist[com[, 1], com[, 1]], 2, function(x) limit.sim(x, par.limit))
      } else {
        prob.death <- prob.death + apply(tr.dist[com[, 1], com[, 1]], 2, limit.sim)
      }
      
      if(any(prob.death<0) | max(prob.death)==0) {
        stop("Death probability cannot be negative, and the maximum must be 
               strictly positive")
      }
    }
    
    # Influence of habitat filtering on mortality
    if("death" %in% type.filt & !is.null(filt)) { 
      
      com_filter <- unlist(env_filter(com[, -(1:2), drop = FALSE]))
      com_filter <- com_filter/max(com_filter)
      
      if (any(is.na(com_filter))) {
        stop("NA values in habitat filter", call. = FALSE)
      }
      
      prob.death <- prob.death + (1 - com_filter)
    }
    
    # If communities contained several traits prob.death object is a one column
    # matrix that needs to be converted to a flat vector for further
    # transformation
    if (is.matrix(prob.death) | is.data.frame(prob.death)) {
      prob.death <- prob.death[, 1]
    }
    
    # Giving names to prob.death
    names(prob.death) <- com[, 1]      
    
    # Position of dead individuals in prob.death vector
    if(d > sum(prob.death > 0)){
      warning("Sampling with replacement in local community (death)", call. = FALSE)
      died <- sample(J, d, replace = TRUE, prob = prob.death)
    } else 
      died <- sample(J, d, replace = FALSE, prob = prob.death)
    
    com <- com[-died, ] 
    
    if (sum(is.na(com[, 1])) != 0) {
      stop("NA values in community composition (2)", call. = FALSE)
    }
  } 
  
  immigrated <- runif(d) < m
  
  # Number of individuals drawn from the regional pool
  J1 <- sum(immigrated)
  # Number of individuals drawn from local offspring
  J2 <- sum(!immigrated)
  
  if (J1 > 0) { 
    # Immigrant drawn from regional pool
    
    # Default equal probability of immigration
    prob = rep(1, nrow(pool))
    
    # Influence of habitat filtering on immigration
    if ("immig" %in% type.filt & !is.null(filt)) {
      pool_filter <- unlist(env_filter(pool[, -(1:2), drop = FALSE]))
      
      if (any(is.na(pool_filter))) {
        stop("NA values in habitat filter", call. = FALSE)
      }
      prob <- prob * pool_filter
    }
    
    # Influence of limiting similarity on immigration
    if("immig" %in% type.limit & !is.null(limit.sim)) {
      # Influence of limiting similarity on establishment
      if(!is.null(limit.sim)) if(length(formals(limit.sim))>1)
      {
        prob.estab <- apply(tr.dist[pool[, 1],
                                   com[, 1]], 1, function(x) limit.sim(x, par.limit))
      } else {
        prob.estab <- apply(tr.dist[pool[, 1],
                                    com[, 1]], 1, limit.sim)
      }
      prob.estab <- prob.estab/max(prob.estab)
      prob.estab <- 1 - prob.estab
      # Add baseline probability
      prob.estab <- prob.estab + 1/J
      prob <- prob * prob.estab
    }
    
    # Add new immigrated individual to community
    if(J1 < sum(prob > 0) & !m.replace) {
      warning("Sampling with replacement from the pool", call. = FALSE)
      com <- rbind(com, pool[sample(1:nrow(pool), J1, replace = TRUE,
                                  prob = prob),])
    } else {
      com <- rbind(com, pool[sample(1:nrow(pool), J1, replace = m.replace,
                                    prob = prob),])
    }
    
    if (any(is.na(com[, 1]))) {
      stop("NA values in community composition (3)", call. = FALSE)
    }
  }
  
  if (J2 > 0) { # Recruitment from com.init
    
    # Default equal probability of recruitment of local offspring
    prob = rep(1, nrow(com.init))
    
    # Influence of habitat filtering on recruitment
    if("loc.recr" %in% type.filt & !is.null(filt)) {
      
      com_filter <- unlist(env_filter(com.init[, -(1:2), drop = FALSE]))
      
      if (any(is.na(com_filter))) {
        stop("NA values in habitat filter", call. = FALSE)
      }
      prob <- prob * com_filter
    }
    
    # Influence of limiting similarity on recruitment
    if("loc.recr" %in% type.limit & !is.null(limit.sim)) {
      # Same constraint is used for local offspring and immigrants
      prob.estab <- prob.death/max(prob.death)
      prob.estab <- 1 - prob.estab
      # Add baseline probability - to be discussed/improved
      prob.estab <- prob.estab + 1/J
      prob <- prob * prob.estab
    }
    
    # Add new established offspring individual to community
    if(J2 < sum(prob < 0)) {
      warning("Sampling with replacement in local community (recruitment)", call. = FALSE)
      com <- rbind(com, com.init[sample(1:nrow(com.init), J2, replace = TRUE,
                                        prob = prob),])
    } else {
      com <- rbind(com, com.init[sample(1:nrow(com.init), J2, replace = FALSE,
                                        prob = prob),])
    }
        
    if (any(is.na(com[, 1]))) {
      print(J2)
      stop("NA values in community composition (4)", call. = FALSE)
    }
  }
  
  if (!is.null(limit.sim)) {
    # If there limiting similarity return the factor
    return(list(com = com, dist.t = dist.t))
    
  } else {
    # Without limiting similarity
    return(list(com = com))
    
  }
}

gauss_limit <- function(dist, par)
{
  sum(exp(-dist^2/(2*par^2)), na.rm = T)
}

