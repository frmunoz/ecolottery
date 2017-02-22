# Function to compute forward simulation of community dynamics with (eventually)
# environmental filtering
forward <- function(initial, prob = 0, d = 1, gens = 150, keep = FALSE,
                    pool = NULL, limit.sim = F, coeff.lim.sim = 1, sigm = 0.1,
                    filt = NULL, prob.death = NULL, method.dist = "euclidean", plot_gens = FALSE) {
  # The function will stop if niche - based dynamics is requested, but trait
  # information is missing in the local community
  # For strictly neutral communities, a vector of species names is enough for
  # the initial community
  
  # Checking basic parameters

  if (!is.numeric(prob) | prob < 0){
    stop("Probability of migration or mutation must be a number belonging to [0; 1] interval.")
  }
  
  if (!is.numeric(d) | d < 0){
    stop("Number of individuals that die in each time step must be a positive number.")
  }
  
  if (!is.numeric(gens) | gens <= 0){
    stop("Number of generations must be a positive number.")
  }
  
  if (!is.logical(keep)){
    stop("keep parameter must be a boolean.")
  }
  
  if (!is.logical(limit.sim)){
    stop("limiting similarity parameter must be a boolean.")
  }
  
  if (!is.numeric(coeff.lim.sim) | coeff.lim.sim < 1){
    stop("coeff.lim.sim parameter must be a positive integer superior to 0.")
  }
  
  if (!is.numeric(sigm) | sigm < 0){
    stop("sigm parameter must be a positive number.")
  }
  
  if (!is.null(filt)){
    if (!is.function(filt)){
      stop("filt must be a function.")
    }
  }    
  
  if ((method.dist %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))==FALSE){
    stop("Provided distance is not existing. See stats::dist function for help.")
  }
  
  if (!is.logical(plot_gens)){
    stop("plot_gens parameter must be a boolean.")
  }
  
  if (plot_gens & !keep){
    warning("plot_gens is valid only if keep parameter is set as TRUE.")
  }
  
  # Stops if only a vector of species name is given as initial community with
  # environmental filtering or limiting similarity
  if ((is.character(initial) | is.vector(initial)) & (limit.sim | !is.null(filt))) {
    stop(paste0("Trait information must be provided along with species",
                " identity in the initial community for niche - based dynamics"))
  }
  
  # If environmental filtering or limiting similarity, the initial community
  # needs to be a matrix or a data.frame
  if (!is.matrix(initial) & !is.data.frame(initial) & (limit.sim | !is.null(filt))) {
    stop("Misdefined initial community")
  }
  
  # If no limiting similarity nor environmental filter -> community dynamics are
  # considered neutral
  if (!limit.sim & is.null(filt)) {
    cat("Simulation of a neutral community\n")
  }
  
  # "pool" will be a three - column matrix of individuals in the regional pool,
  # with individual id in first column, species name in second column, and
  # additional trait information for niche - based dynamics in third column
  if (is.character(pool)) {
    pool <- data.frame(id = 1:length(pool),
                       sp = pool,
                       trait = rep(NA, length(pool)),
                       stringsAsFactors = F)
  
    if (limit.sim | !is.null(filt)) {
      cat("No trait information provided in the regional pool\n")
      pool[, 3] <- runif(nrow(pool))
      cat("Random trait values attributed to individuals of the regional pool\n")
      colnames(pool) <- c("id", "sp", "trait")
    }
  }
  
  # If species pool is specified by user
  if (!is.null(pool)) {
    if (ncol(pool) < 2) {
      stop(paste0("The regional pool is misdefined (at least two columns ",
                  "required when a matrix or data frame is provided)"))
    } else if (ncol(pool) == 2) {
      cat("No trait information provided in the regional pool\n")
    }
    if ((!limit.sim | !is.null(filt)) & ncol(pool) < 3) {
      pool[, 3] <- runif(nrow(pool))
      
      cat(paste0("Random (uniform) trait values attributed to individuals of ",
                 "the regional pool\n"))
    }
      colnames(pool) <- c("id", "sp", "trait")
  }
  
  # "init_comm" is a 3 columns matrix of individuals in the initial community,
  # with individual id in first column, species name in second column, and
  # additional trait information for niche-based dynamics in the third column
  if (is.character(initial)) {  # If only list of species names provided
    J <- length(initial)
    # The ids of individuals present in the initial community begi with "init"
  init_comm <- data.frame(id = paste("init", 1:J, sep = ""),
                          sp = initial,
                          trait = rep(NA, J),
                          stringsAsFactors = F) 
  } else {
    
    if (ncol(initial) < 3) {
      cat(paste0("Two-column initial community: assumed to represent species ",
                 "and trait information; individual ids will be generated"))
    J <- nrow(initial)
    init_comm <- data.frame(id = paste("init", 1:J, sep = ""),
                            sp = initial[, 1],
                            trait = initial[, 2],
                            stringsAsFactors = F)
	    } else
	    {
	    	init_comm <- initial
	    }
  }
  
  if ((limit.sim | !is.null(filt)) & any(is.na(init_comm[, 3]))) {
    stop(paste0("Trait information must be provided in initial community",
                "composition for niche-based dynamics"))
  }
  
  colnames(init_comm) <- c("id", "sp", "trait")
  
  new.index <- 0
   
  ## Forward simulation with community
  
  # Begins with the initial community
  next_comm <- init_comm
  
  if (keep) {  # If the user asked to keep all the communities at each timestep
    comm_through_time <- c()
    limit.sim.t <- c()
    
    # Simulate the community for the given number of generations
    for (i in 1:gens) {
      comm_through_time[[i]] <- next_comm  # Store the community at time i
      
      # Simulate community dynamics
      next_comm <- pick(next_comm, d = d, prob = prob, pool = pool,
                        prob.death = prob.death, limit.sim = limit.sim,
                        coeff.lim.sim = coeff.lim.sim, sigm = sigm,
                        filt = filt, new.index = new.index, method.dist = "euclidean")
      
      # Store limiting similarity matrix if simulated
      if (!is.null(limit.sim)) {
        limit.sim.t <- c(limit.sim.t, next_comm$limit.sim.t)
      }
      new.index <- next_comm$new.index
      next_comm <- next_comm$com
    } 
    
    if (plot_gens) { # Plotting number of individuals and species over generations
      
      uniq_list <- lapply(comm_through_time,
                          function(y) apply(y, 2, function(x) length(unique(x))))
      
      uniq_df <- do.call(rbind.data.frame, uniq_list)
      colnames(uniq_df) <- c("id_uniq", "nb_sp", "nb_tra")
      uniq_df$gens <- seq(1:nrow(uniq_df))
      uniq_df$nb_id <- J
      
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        
        # Plot the number of individuals through all the generations
        plot_individuals = ggplot2::ggplot(uniq_df, ggplot2::aes(gens, nb_id)) +
                ggplot2::geom_line() +
                ggplot2::geom_line(aes(gens, id_uniq), size = 1) +
                ggplot2::labs(x = "Number of generations",
                              y = "Number of unique individuals")
        
        # Plot the number of species through all the generations
        plot_species = ggplot2::ggplot(uniq_df, ggplot2::aes(gens, nb_sp)) +
                ggplot2::geom_line(size = 1) +
                ggplot2::labs(x = "Number of generations",
                     y = "Number of unique species")
        
        print(plot_individuals)
        print(plot_species)
      }
      
    }
    
    return(list(com_t = comm_through_time,
                limit.sim.t = limit.sim.t,
                pool = pool))
    
  } else {# Keep only the last community
    for (i in 1:gens) {
      # Simulate community dynamics for a timestep
      next_comm <- pick(next_comm, d = d, prob = prob, pool = pool,
                        prob.death = prob.death, limit.sim = limit.sim,
                        coeff.lim.sim = coeff.lim.sim, sigm = sigm,
                        filt = filt, new.index = new.index, method.dist = "euclidean")
      
      new.index <- next_comm$new.index
      next_comm <- next_comm$com
    }
    return(list(com = next_comm, pool = pool))
  }
}

# Precise function to simulate a single timestep by picking an individual in
# the pool or make an individual mutate
pick <- function(com, d = 1, prob = 0, pool = NULL, prob.death = prob.death,
                 limit.sim = NULL, coeff.lim.sim = 1, sigm = 0.1, filt = NULL,
                 new.index = new.index, method.dist = "euclidean") {
  
  
  if (is.null(pool)) {
    # If no species pool specified, mutate an individual
    return(pick.mutate(com, d = d, prob.of.mutate = prob, new.index = new.index)) 
  
  } else {
	  
    if((!is.null(filt) | limit.sim) & prob > 0 & any(is.na(pool[,3]))) {
	    stop("With environmental filtering, NA trait values not allowed in regional pool")
    }
	  
  # If there is a species pool make an individual immigrates
    return(pick.immigrate(com, d = d, prob.of.immigrate = prob, pool = pool,
                          prob.death = prob.death, limit.sim = limit.sim,
                          coeff.lim.sim = coeff.lim.sim, sigm = sigm,
                          filt = filt, method.dist = "euclidean"))
  }
}

# Return community with mutated inidividual (= new species)
pick.mutate <- function(com, d = 1, prob.of.mutate = 0, new.index = 0) {
  
  if (is.vector(com)) {
    # If community only defined by species names
    J <- length(com)
    
    com <- data.frame(id = paste("ind", 1:J, sep = ""),
                      sp = as.character(com),
                      trait = rep(NA, J),
                      stringsAsFactors = F)
  
  } else if (is.matrix(com) | is.data.frame(com)) {
    # If the community has defined traits
    J <- nrow(com)
    
    if (is.matrix(com)) {
      	com <- as.data.frame(com, stringsAsFactors = F)
 	com[, 1] <- as.character(com[, 1])
     	com[, 2] <- as.character(com[, 2])	  
    }
    
  } else {
    stop("pick.mutate: misdefined community composition")
  }
  
  ## Simulate the dynamics of the community
  
  # Number of individuals who die at this timestep
  died <- sample(J, d, replace = TRUE)
  
  # How many of the dead individuals are replaced by mutated individuals
  mutated <- runif(length(died)) < prob.of.mutate
  
  # Number of mutated individuals
  n_mutated <- sum(mutated)
  
  # Dead individuals who did not mutate
  dead_non_mutated <- sum(!mutated)
  
  # Replace dead non mutated individuals by individuals from other species
  com[died[!mutated], ] <- com[sample(1:nrow(com), dead_non_mutated,
                                      replace = TRUE), ]
  
  # When mutation occurs
  if (n_mutated > 0) {
    
    # Attribute new species to individuals who mutate
    com[died[mutated], 1:2] <- paste("new.sp", new.index + (1:n_mutated),
                                     sep = "")
    
    #com[died[mutated], 3] <- rep(NA, n_mutated)  # No trait values
    # Default to be a trait value drawn from uniform distribution between 0 and 1
    com[died[mutated], 3] <- runif(n_mutated)
	  
    # Number of new species which appeared (next one will be new.index + 1)
    new.index <- new.index + n_mutated
  }
  
  return(list(com = com, new.index = new.index))
}

# Function to return individuals who immigrated from the species pool
# limit.sim = distances de traits; filt = habitat filtering function
pick.immigrate <- function(com, d = 1, prob.of.immigrate = 0, pool,
                           prob.death = NULL, limit.sim = NULL,
                           coeff.lim.sim = 1, sigm = 0.1, filt = NULL, method.dist = "euclidean") {
  
  if (is.vector(com)) {
    # If community only defined by species names
    J <- length(com)
    
    com <- data.frame(id = paste("ind", 1:J, sep = ""),
                    sp = as.character(com),
                    trait = rep(NA, J),
                    stringsAsFactors = F)
  
  } else if (is.matrix(com) | is.data.frame(com)) {
    # If the community has defined traits
    J <- nrow(com)
    
    if (is.matrix(com)) {
      com <- as.data.frame(com, stringsAsFactors = F)
    }
    com[, 1] <- as.character(com[, 1])
    com[, 2] <- as.character(com[, 2])
      
  } else {
    stop("pick.immigrate: misdefined community composition")
  }
  
  # Function defining habitat filtering according to trait value
  if (!is.null(filt)) {
    hab_filter <- function(x) filt(x)
  } else {
    # If no function defined, dummy function returning one
	  hab_filter <- function(x) sapply(x, function(x) 1)
  }
  
  # Limiting similarity depends on community composition at each time step
  if (limit.sim) {
   limit.sim <- as.matrix(dist(com[, 3], method = method.dist))
   colnames(limit.sim) <- com[, 1]
   rownames(limit.sim) <- com[, 1]
   diag(limit.sim) <- NA
  } else {
    limit.sim <- NULL
  }
					   
  # Initial community
  com.init <- com
  
  if (is.null(limit.sim) & is.null(filt)) {
    
    died <- sample(J, d, replace = TRUE)
    
    com <- com[-died, ]
    
    if (any(is.na(com[, 1]))) {
      stop("Error: NA values in community composition (1)")
      }
  
  } else {# Influence of limiting similarity and habitat filtering on mortality
    
    # Vector of the individual probability of dying
    if (is.null(prob.death)) {
      prob.death <- rep(1, nrow(com))
    }
    
    # Under limiting similarity, mortality increases when an individual is more
    # similar to other resident individuals
    if (!is.null(limit.sim)) {
      
      if (sum(!com[, 1] %in% rownames(limit.sim)) > 0 |
          sum(!com[, 1] %in% colnames(limit.sim)) > 0) {
        stop("limit.sim: mismatch of species names")
      }
      
      # For each species: compute limiting similarity coefficient based on Gaussian distribution
      limit.sim.t <- apply(limit.sim[com[, 1], com[, 1]], 2,
                           function(x) (sum(exp( -x^2 / (2*(sigm^2))), na.rm = T)))
      
      prob.death <- coeff.lim.sim*limit.sim.t
      # Scaling prob.death
      prob.death <- (prob.death - min(prob.death)) / (max(prob.death) - min(prob.death))

      # If all probabilities null, sample won't work. An identical and weak probability is given to each species.
      if(sum(prob.death)==0){
        prob.death <- prob.death + 0.001
      }
      
      # limit.sim.t will display the average distance between trait of species for the whole community
      limit.sim.t <- mean(limit.sim[com[, 1], com[, 1]], na.rm = T)
      
    }
    # Habitat filtering also influences the individual death probability
    if (!is.null(filt)) {
      if (any(is.na(hab_filter(com[, 3])))) {
        stop("Error: NA values in habitat filter")
      }
      prob.death <- prob.death * (1 - hab_filter(com[, 3]) / sum(hab_filter(com[, 3])))
      
      # Giving names to prob.death
      names(prob.death) <- com[, 1]
      
      # If all probabilities null, sample won't work. An identical and weak probability is given to each species.
      if(sum(prob.death)==0){
        prob.death <- prob.death + 0.001
      }
    }
    
    # Position of dead individuals in prob.death vector
    died <- sample(J, d, replace = T, prob = prob.death)
    # Identities of dead individuals
    id_died <- names(prob.death)[died]
    # If several individuals concerned: the first one of each identity dies
    id_died_pos <- as.numeric(sapply(id_died, function(x) which(com[, 1] %in% x)[1]))
    # Community composition after mortality: died individuals removed
    com <- com[-id_died_pos, ] 
    
    if (sum(is.na(com[, 1])) != 0) {
      stop("Error: NA values in community composition (2)")
    }
  } 
  
  immigrated <- runif(d) < prob.of.immigrate
  
  # If probability of immigration is high, then the new individual is drawn from the regional pool
  J1 <- sum(immigrated)
  # The lower the probability of immigration, the higher the probability of drawing the new individual from the community
  J2 <- sum(!immigrated)
  
  if (J1 > 0) { # Immigrant drawn from regional pool
    # Influence of habitat filtering on immigration
    if (is.null(filt)) {
      add <- pool[sample(1:nrow(pool), J1, replace = TRUE), 1:3]
      com <- rbind(com, add)
      
      if (sum(is.na(com[, 1])) != 0) {
        stop("Error: NA values in immigrants")
      }
    } else {
      com <- rbind(com, pool[sample(1:nrow(pool), J1, replace = TRUE,
                                prob = hab_filter(pool[, 3])), 1:3])
      
      if (any(is.na(hab_filter(pool[, 3])))) {
        stop("Error: NA values in habitat filtering")
      }
    }
    
    if (any(is.na(com[, 1]))) {
      stop("Error: NA values in community composition (3)")
      }
  }
  
  if (J2 > 0) { # Immigrant drawn from com.init
    
    com <- rbind(com, com.init[sample(1:nrow(com.init), J2, replace = TRUE), 1:3])
    
    if (any(is.na(com[, 1]))) {
      print(J2)
      stop("Error: NA values in community composition (4)")
    
    }
  }
  
  
  if (!is.null(limit.sim)) {
    # If there limiting similarity return the factor
    return(list(com = com, limit.sim.t = limit.sim.t))
  
  } else {
    # Without limiting similarity
    return(list(com = com))
  
  }
}


