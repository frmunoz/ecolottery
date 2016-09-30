forward <- function(initial, prob = 0, D = 1, gens = 150, keep = FALSE,
                    pool = NULL, limit.sim = F, coeff.lim.sim = 2, sigma = 0.1,
                    filt = NULL, prob.death = NULL, method.dist = "euclidean") {
  # The function will stop if niche - based dynamics is requested, but trait information is missing in the local community
  # For strictly neutral communities, a vector of species names is enough for the initial community
  if ((is.character(initial) | is.vector(initial)) & (limit.sim | !is.null(filt))) {
    stop("Trait information must be provided along with species identity in the initial community for niche - based dynamics")
  }
  
  if (!is.matrix(initial) & !is.data.frame(initial) & (limit.sim | !is.null(filt))) {
    stop("Misdefined initial community")
  }
  
  if (!limit.sim & is.null(filt)) {
    cat("Simulation of a neutral community\n")
  }
  
  # "pool" will be a three - column matrix of individuals in the regional pool, with individual id in first column, species name in second column, and additional trait information for niche - based dynamics in third column
  if (is.character(pool)) {
    pool <- data.frame(id = 1:length(pool), sp = pool, trait = rep(NA, length(pool)), stringsAsFactors = F)
  
    if (limit.sim | !is.null(filt)) {
      cat("No trait information provided in the regional pool\n")
      pool[, 3] <- runif(nrow(pool))
      cat("Random trait values attributed to individuals of the regional pool\n")
      colnames(pool) <- c("id", "sp", "trait")
    }
  }
  
  if (!is.null(pool)) {
    if (ncol(pool)<2) {
      stop("The regional pool is misdefined (at least two columns required when a matrix or data frame is provided)")
    } else if (ncol(pool) == 2) {
      cat("No trait information provided in the regional pool\n")
    }
    if (!limit.sim | !is.null(filt)) {
      pool[, 3] <- runif(nrow(pool))
      cat("Random trait values attributed to individuals of the regional pool\n")
      colnames(pool) <- c("id", "sp", "trait")
    }
  }
  
  # "a" will be a three - column matrix of individuals in the initial community, with individual id in first column, species name in second column, and additional trait information for niche - based dynamics in third column
  if (is.character(initial)) {
    J <- length(initial)
    # The ids of individuals present in the initial community begi with "init"
  a <- data.frame(id = paste("init", 1:J, sep = ""), sp = initial, trait = rep(NA, J), stringsAsFactors = F) 
  } else {
    if (ncol(initial) < 3) {
      cat("Two - column initial community: assumed to represent species and trait information; individual ids will be generated")
    }
  	J <- nrow(initial)
    a <- data.frame(id = paste("init", 1:J, sep = ""), sp = initial[, 1], trait = initial[, 2], stringsAsFactors = F)
  }
  
  if ((limit.sim | !is.null(filt)) & any(is.na(a[, 3]))) {
    stop("Trait information must be provided in initial community composition for niche - based dynamics")
  }
  colnames(a) <- c("id", "sp", "trait")
  
  #if(limit.sim | !is.null(filt)) 
  # {
  # sp.pool.tmean <- tapply(pool[, 3], pool[, 2], mean)
  # sp.init.tmean <- tapply(a[, 3], a[, 2], mean)
  # if(any(sp.pool.tmean[rownames(sp.init.tmean)]!= sp.init.tmean) | any(sp.init.tmean[rownames(sp.pool.tmean)]!= sp.pool.tmean)) warning("Mismatch of mean trait values per species in initial community and in regional pool")
  # }
  
  # Limiting similarity will be based on trait distances in the metacommunity plus the initial community
  if (limit.sim) {
    if (!is.null(pool)) {
      limit.sim <- as.matrix(dist(c(a[, 3], pool[, 3]), method = method.dist))
      colnames(limit.sim) <- c(a[, 1], pool[, 1])
      rownames(limit.sim) <- c(a[, 1], pool[, 1])
    } else {
      limit.sim <- as.matrix(dist(a[, 3], method = method.dist))
      colnames(limit.sim) <- a[, 1]
      rownames(limit.sim) <- a[, 1]
    }
    diag(limit.sim) <- 0
  } else {
    limit.sim <- NULL
  }
  
  new.index <- 0
  
  if (keep) {
    aa <- c()
    limit.sim.t <- c()
    for (i in 1:gens) {
      aa[[i]] <- a
      a <- pick(a, D = D, prob = prob, pool = pool, prob.death = prob.death, limit.sim = limit.sim, coeff.lim.sim = coeff.lim.sim, sigma = sigma, filt = filt, new.index = new.index)
      if (!is.null(limit.sim)) {
        limit.sim.t <- c(limit.sim.t, a$limit.sim.t)
      }
      new.index <- a$new.index
      a <- a$a
    } 
    return(list(aa = aa, limit.sim.t = limit.sim.t, pool = pool))
  } else {
    for (i in 1:gens) {
      a <- pick(a, D = D, prob = prob, pool = pool, prob.death = prob.death,
                limit.sim = limit.sim, coeff.lim.sim = coeff.lim.sim,
                sigma = sigma, filt = filt, new.index = new.index)
      new.index <- a$new.index
      a <- a$a
    }
    return(list(a = a, pool = pool))
  }
}
  
# Return "mutate" individuals if pool = NULL, else return immigrate individuals
pick <- function(a, D = 1, prob = 0, pool = NULL, prob.death = prob.death,
                 limit.sim = NULL, coeff.lim.sim = 2, sigma = 0.1, filt = NULL,
                 new.index = new.index) {
  if(is.null(pool)) {
    return(pick.mutate(a, D = D, prob.of.mutate = prob, new.index = new.index)) 
  } else {
	return(pick.immigrate(a, D = D, prob.of.immigrate = prob, pool = pool,
	                      prob.death = prob.death, limit.sim = limit.sim,
	                      coeff.lim.sim = coeff.lim.sim, sigma = sigma,
	                      filt = filt))
  }
}

# Return "mutate" individuals
pick.mutate <- function(a, D = 1, prob.of.mutate = 0, new.index = 0) {
  if (is.vector(a)) {
    J <- length(a)
    a <- data.frame(id = paste("ind", 1:J, sep = ""), sp = as.character(a), trait = rep(NA, J), stringsAsFactors = F)
  }
  else if (is.matrix(a) | is.data.frame(a)) {
    J <- nrow(a)
    if (is.matrix(a)) {
      a <- as.data.frame(a, stringsAsFactors = F)
      a[, 1] <- as.character(a[, 2])
    }
  } else {
    stop("pick.mutate: misdefined community composition")
  }
  
  died <- sample(J, D, replace = TRUE)
  mutated <- runif(length(died)) < prob.of.mutate
  J1 <- sum(mutated)
  J2 <- sum(!mutated)
  a[died[!mutated], ] <- a[sample(1:nrow(a), J2, replace = TRUE), ]
  if (J1 > 0) {
    # When a mutation occurs, the same individual id and species name are provided
    a[died[mutated], 1:2] <- paste("new.sp", new.index + (1:J1), sep = "")
    a[died[mutated], 3] <- rep(NA, J1)
    new.index <- new.index + J1
  }
  return(list(com = a, new.index = new.index))
}

# Return "immigrate" individuals
# limit.sim = distances de traits; filt = fonction pour habitat filtering
pick.immigrate <- function(a, D = 1, prob.of.immigrate = 0, pool, prob.death = NULL, limit.sim = NULL, coeff.lim.sim = 2, sigma = 0.1, filt = NULL) {
  if (is.vector(a)) {
    J <- length(a)
    a <- data.frame(id = paste("ind", 1:J, sep = ""), sp = as.character(a), trait = rep(NA, J), stringsAsFactors = F)
  }
  else if (is.matrix(a) | is.data.frame(a)) {
    J <- nrow(a)
    if (is.matrix(a)) {
      a <- as.data.frame(a, stringsAsFactors = F)
      a[, 2] <- as.character(a[, 2])
    }
  } 
  else stop("pick.immigrate: misdefined community composition")
  
  # Function defining habitat filtering according to trait value
  if (!is.null(filt)) {
    filt1 <- function(x) filt(x)
  } else {
	  filt1 <- function(x) sapply(x, function(x) 1)
  }
  
  a.init <- a
  # For debug:
  #print(a[(nrow(a) - 10):nrow(a), ])
  #print(head(pool))
  if (is.null(limit.sim) & is.null(filt)) {
    died <- sample(J, D, replace = TRUE)
    a <- a[ - died, ]
    if (any(is.na(a[, 1]))) {
      stop("Error: NA values in community composition (1)")
      }
  } else {# Influence of limiting similarity and habitat filtering on mortality
    # Vector of the individual probability of dying
    if (is.null(prob.death)) {
      prob.death <- rep(1, nrow(a))
    }
    # Under limiting similarity, mortality increases when an individual is more similar to other resident individuals
    if (!is.null(limit.sim)) {
      if (sum(!a[, 1] %in% rownames(limit.sim)) > 0 | sum(!a[, 1] %in% colnames(limit.sim)) > 0) {
        stop("limit.sim: mismatch of species names")
      }
      limit.sim.t <- apply(limit.sim[a[, 1], a[, 1]], 2, function(x) (sum(exp( - x^2 / (2*sigma^2)))))
      prob.death <- (coeff.lim.sim - 1)*limit.sim.t
      limit.sim.t <- mean(limit.sim.t)
    }
    # Habitat filtering also influences the individual death probability
    if (!is.null(filt)) {
      if (any(is.na(filt1(a[, 3])))) {
        stop("Error: NA values in habitat filtering")
      }
      prob.death <- prob.death * (1 - filt1(a[, 3]) / sum(filt1(a[, 3])))
    }
    died <- rmultinom(1, D, prob.death)
    a <- a[ -died, ] # Community composition after mortality
    if (sum(is.na(a[, 1])) != 0) {
      stop("Error: NA values in community composition (2)")
    }
  } 
  immigrated <- runif(D) < prob.of.immigrate
  J1 <- sum(immigrated)
  J2 <- sum(!immigrated)
  if (J1 > 0) {
    # Influence of habitat filtering on immigration
    if (is.null(filt)) {
      add <- pool[sample(1:nrow(pool), J1, replace = TRUE), 1:3]
      a <- rbind(a, add)
      if (sum(is.na(a[, 1])) != 0) {
        stop("Error: NA values in immigrants")
      }
    } else {
      a <- rbind(a, pool[sample(1:nrow(pool), J1, replace = TRUE, prob = filt1(pool[, 3])), 1:3])
      if (any(is.na(filt1(pool[, 3])))) {
        stop("Error: NA values in habitat filtering")
      }
    }
    if (any(is.na(a[, 1]))) {
      stop("Error: NA values in community composition (3)")
      }
  }
  if (J2 > 0) {
    a <- rbind(a, a.init[sample(1:nrow(a.init), J2, replace = TRUE), 1:3])
    if (any(is.na(a[, 1]))) {
      print(J2)
      stop("Error: NA values in community composition (4)")
    }
  }
  if (!is.null(limit.sim)) {
    return(list(a = a, limit.sim.t = limit.sim.t))
  } else {
    return(list(a = a))
  }
}


