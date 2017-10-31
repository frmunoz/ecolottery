coalesc_abc <- function(comm.obs, pool = NULL, multi = "single", prop = F, traits = NULL,
                        f.sumstats, filt.abc = NULL, params = NULL,
                        theta.max = NULL, nb.samp = 10^6, parallel = TRUE,
                        tol = NULL, pkg = NULL, method = "rejection")
{
  
  if(!method%in%c("rejection", "loclinear", "neuralnet", "ridge"))
    stop("method.abc should be either rejection, loclinear, neuralnet or ridge")
  
  if (!is.function(f.sumstats)) {
    stop("You must provide a function to calculate summary statistics",
         "(f.sumstats)")
  }
  
  if(prop & multi!="tab")
    stop("prop data can only be handled in tab format")
  
  if (length(formals(f.sumstats)) > 2) {
    stop("f.sumstats must be a function of up to two arguments")
  }
  
  if (is.null(tol)){
    warning("You must provide a tolerance value for ABC computation.\n",
            "The function will only provide simulations and will not perform ",
            "ABC analysis.")
  } else if(is.null(params)){
    warning("No value provided for params argument. Only m and theta will be ",
            "estimated.")
  }
  
  if (is.character(comm.obs)) {
      comm.obs <- data.frame(id = 1:length(comm.obs),
                       sp = comm.obs,
                       trait = rep(NA, length(comm.obs)),
                       stringsAsFactors = FALSE)
      colnames(comm.obs) <- c("id", "sp", "trait")
  }
  
  if (is.null(pool)) {
    warning(paste0("No species pool provided: pool will be simulated ",
                   "and logseries theta estimated"))
    if (is.null(theta.max))
      theta.max <- 500
  }
  
  if (!requireNamespace("abc", quietly = TRUE)) {
    stop("coalesc_abc requires package abc to be installed")
  }
  
  # Other required packages
  for (i in pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package ", pkg, " is not available", sep = ""))
    }   
  }
  
  # Community size
  if (!(multi %in% c("single", "tab", "seqcom"))){
    stop("multi parameter must be either single, tab or seqcom.")
  }
  
  if (multi == "single") {
    J <- nrow(comm.obs)
    nb.com <- 1
  } else {
    if (multi == "tab") {
      # comm.obs is a species-by-site matrix/data.frame
      J <- apply(comm.obs, 1, function(x) sum(x, na.rm = TRUE))
      # if the dataset includes relative proportions, the columns must sum to 1
      if(prop & any(J!=1)) stop("Relative species abundances must sum to 1")
      nb.com <- nrow(comm.obs)
    } else if (multi == "seqcom") {
      # comm.obs is a list of communities with individuals on rows in each
      # community
      J <- lapply(comm.obs, nrow)
      nb.com <- length(comm.obs)
    }
  }
  
  # Trait values can be provided with community composition
  # Mean trait values of pool are stored in traits in absence of trait
  # information in local community
  if (!is.null(pool)){
    if(ncol(pool) >= 3) {
      traits <- data.frame(apply(data.frame(pool[,-(1:2)]), 2,
                                 function(x) {
                                   tapply(x, pool[, 2],
                                          function(y)
                                            mean(y, na.rm = TRUE)
                                   )}))
    }
  }
  
  if (is.null(traits)){
    warning("Trait information is not provided")
  }
  
  if(multi == "seqcom"){
    if (length(formals(f.sumstats))==1) {
      stats.obs <- lapply(comm.obs, f.sumstats)
    } else {
      stats.obs <- lapply(comm.obs, function(x) f.sumstats(x, traits))
    }
  } else {
    if (length(formals(f.sumstats))==1) {
      stats.obs <- f.sumstats(comm.obs)
    } else {
      stats.obs <- f.sumstats(comm.obs, traits)
    }
  }
  
  # Community simulation
  sim <- do.simul(J, pool, multi, prop, nb.com, traits, f.sumstats, filt.abc, params,
                  theta.max, nb.samp, parallel, tol, pkg, method)
  
  # Scaling f.sumstats criterions for simulations
  if (multi == "seqcom"){
    stats.mean <- lapply(sim$stats, function(x) {
      apply(x, 2, function(y) mean(y, na.rm = TRUE))
    })
    
    stats.sd <- lapply(sim$stats, function(x) {
      apply(x, 2, function(y) sd(y, na.rm = TRUE))
    })
    
    sim$stats.scaled <- list()
    for (i in 1:length(sim$stats)) {
      sim$stats.scaled[[i]] <- t(apply(sim$stats[[i]], 1, function(x) {
        (x - stats.mean[[i]])/stats.sd[[i]]
      }))
      colnames(sim$stats.scaled[[i]]) <- colnames(sim$stats[[i]])
    }
  } else {
    stats.mean <- apply(sim$stats, 2, function(x) mean(x, na.rm = TRUE))
    stats.sd <- apply(sim$stats, 2, function(x) sd(x, na.rm = TRUE))
    sim$stats.scaled <- t(apply(sim$stats, 1,
                                function(x) (x - stats.mean)/stats.sd))
    colnames(sim$stats.scaled) <- colnames(sim$stats)
  }
  
  if (is.null(tol)){
    res.abc <- NA
  } else {
    if (multi == "seqcom") {
      stats.obs.scaled <- list()
      for (i in 1:length(sim$stats)){
        stats.obs.scaled[[i]] <- (stats.obs[[i]] - stats.mean[[i]])/
          stats.sd[[i]]
      }
    } else {
      stats.obs.scaled <- (stats.obs - stats.mean)/stats.sd
    }
    
    if (is.null(tol)){
      res.abc <- NA
    } else {
      if (multi == "seqcom"){
        # ABC estimation
        sel <- lapply(sim$stats.scaled, function(x) {
          which(rowSums(is.na(x)) == 0)
        })
        
        res.abc <- list()
        for(i in 1:length(stats.obs.scaled)){
          res.abc[[i]] <- tryCatch(
            abc::abc(target = stats.obs.scaled[[i]],
                     param = sim$params.sim[[i]][sel[[i]],],
                     sumstat = sim$stats.scaled[[i]][sel[[i]],],
                     tol = tol,
                     method = method),
            error = function(x) {
              warning("ABC computation failed with the requested method.")
            }
        )}
          for(i in 1:length(res.abc)){
            if (is.character(res.abc[[i]])){
              res.abc <- NA
            }
          }
      } else {
        # ABC estimation
        sel <- which(rowSums(is.na(sim$stats.scaled)) == 0)
        
        res.abc <- tryCatch(
          abc::abc(target = stats.obs.scaled,
                   param = sim$params.sim[sel,],
                   sumstat = sim$stats.scaled[sel,],
                   tol = tol,
                   method = method),
          error = function(x) {
            warning("ABC computation failed with the requested method.")
          }
        )
        if (is.character(res.abc)){
          res.abc <- NA
        }
      }
    }
  }
  
  return(list(par = sim$params.sim, obs = stats.obs,
              obs.scaled = stats.obs.scaled, ss = sim$stats.scaled,
              ss.scale = data.frame(mean=stats.mean,sd=stats.sd),
              abc = res.abc))
}

do.simul <- function(J, pool = NULL, multi = "single", prop = F, nb.com = NULL,
                     traits = NULL, f.sumstats = NULL, filt.abc = NULL, params,
                     theta.max = NULL, nb.samp = 10^6, parallel = TRUE,
                     tol = NULL, pkg = NULL, method = "rejection") {
  
  if (!requireNamespace("parallel", quietly = TRUE) & parallel) {
    warning("parallel = TRUE requires package 'parallel' to be installed\n",
            "changed to parallel = FALSE")
    parallel <- FALSE
  }
  
  if (parallel) {
    # Start up a parallel cluster
    parCluster <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
  }
  
  # Uniform prior distributions of parameters
  prior <- c()
  if(!is.null(params)) {
    for (i in 1:nrow(params)) {
      prior[[i]] <- runif(nb.samp, min = params[i, 1], max = params[i, 2])
    }
  }
  names(prior) <- rownames(params)
  
  # Note - Defining a lower bound at 0 for m and theta can entail issues when species richness is 1
  # in simulated community; we should allow the user to define the prior for m and theta in
  # the future
  
  prior[[length(prior) + 1]] <- runif(nb.samp, min = 0, max = 1)
  names(prior)[length(params) + 1] <- "m"
  
  if(prop) 
  {
    prior[[length(prior)+1]] <- runif(nb.samp, min = 100, max = 1000)
    warning("The prior of community size is uniform between 100 and 1000")
  }
  
  if (is.null(pool)) {
    prior[[length(prior) + 1]] <- runif(nb.samp, min = 0, max = theta.max)
    names(prior)[length(prior)] <- "theta"
  }
  
  # Function to perform simulations
  mkWorker <- function(traits, nb.com, multi, prop, prior, J, pool, filt.abc,
                       f.sumstats, pkg) {
    force(J)
    force(pool)
    force(traits)
    force(filt.abc)
    force(f.sumstats)
    force(prop)
    force(prior)
    force(nb.com)
    force(pkg)
    force(multi)
    
    summCalc <- function(j, multi, traits, nb.com, prior, J, prop, pool, filt.abc,
                         f.sumstats) {
      
      params.samp <- unlist(lapply(prior,function(x) x[j]))
      stats.samp <- NA
      params.samp.all <- params.samp
      names(params.samp.all) <- names(prior)
      
      if(prop) if(!is.null(pool))
      {
        J <- round(params.samp[length(params.samp)])
        params.samp <- params.samp[-length(params.samp)]
      } else if(is.null(pool)) 
      {
        J <- round(params.samp[length(params.samp)-1])
        params.samp <- params.samp[-(length(params.samp)-1)]
      }
          
      if (is.null(pool)) {
        if (multi == "seqcom"){
          pool <- coalesc(mean(unlist(J))*100,
                          theta = params.samp[length(params.samp)])$pool
        } else {
          pool <- coalesc(mean(J)*100,
                          theta = params.samp[length(params.samp)])$pool
        }
        params.samp <- params.samp[-length(params.samp)]
      }
      
      if (!is.null(filt.abc)) {
        filt <- function(x) filt.abc(x, params.samp[-length(params.samp)])
      } else {
        filt <- NULL
      }
      
      if (nb.com > 1) {
        if(multi == "tab") {
          pool.sp <- unique(pool$sp)
          meta.samp <- array(0, c(nb.com, length(pool.sp)))
          colnames(meta.samp) <- pool.sp
          
          for (i in 1:nb.com) {
            try({
              comm.samp <- coalesc(J[i], m = params.samp[length(params.samp)],
                                   filt = filt,
                                   pool = pool, traits = traits)
              tab <- table(comm.samp$com[,2])
              meta.samp[i,names(tab)] <- tab
            })
          }
          
          if (length(formals(f.sumstats))==1) {
            stats.samp <- f.sumstats(meta.samp)
          } else {
            stats.samp <- f.sumstats(meta.samp, traits)
            }
        } else if(multi == "seqcom") {
          seqcom.samp <- list()
          
          for (i in 1:nb.com) {
            try({
              seqcom.samp[[i]] <- coalesc(J[[i]],
                                          m = params.samp[length(params.samp)],
                                          filt = filt,
                                          pool = pool, traits = traits)
            })
          }
          if (length(formals(f.sumstats)) == 1) {
            seqcom.samp.com <- lapply(seqcom.samp, function(l) l[[1]])
            stats.samp <- lapply(seqcom.samp.com, f.sumstats)
          } else {
            seqcom.samp.com <- lapply(seqcom.samp, function(l) l[[1]])
            stats.samp <- lapply(seqcom.samp.com, function(x) {
              f.sumstats(x, traits)
            })
          }
        }
      
      } else { # single community
        comm.samp <- coalesc(J,
                             m = params.samp[length(params.samp)],
                             filt = filt,
                             pool = pool, traits = traits)
        if(prop) comm.samp$com <- t(table(comm.samp$com[,2])/J)
        if (length(formals(f.sumstats))==1) {
          stats.samp <- f.sumstats(comm.samp$com)
        } else {
          stats.samp <- f.sumstats(comm.samp$com, traits)
        }
      }
          
      return(list(sum.stats = stats.samp, param = params.samp.all))
    }
    
    worker <- function(j) {
      if (!requireNamespace("ecolottery", quietly = TRUE)) {
        stop(paste("Package ecolottery is not available", sep = ""))
      }
      # Other required packages
      if (!is.null(pkg)) for (i in 1:length(pkg)) {
        if (!requireNamespace(pkg[i], quietly = TRUE)) {
          stop(paste("Package ", pkg[i], " is not available", sep = ""))
        }
      }
      summCalc(j, multi, traits, nb.com, prior, J, prop, pool, filt.abc, f.sumstats)
    }
    return(worker)
  }
  
  # Calculation of summary statistics over the whole range of parameters
  if (parallel) {
    models <- parallel::parLapply(parCluster, 1:nb.samp,
                                  mkWorker(traits, nb.com, multi, prop, prior, J,
                                           pool, filt.abc, f.sumstats, pkg))
  } else {
    models <- lapply(1:nb.samp, mkWorker(traits, nb.com, multi, prop, prior, J, pool,
                                         filt.abc, f.sumstats, pkg))
  }
  
  if (parallel) {
    # Close parCluster
    if (!is.null(parCluster)) {
      parallel::stopCluster(parCluster)
      parCluster <- c()
    }
  }
  
  if (multi == "seqcom"){ # Results stored in a list
    stats <- list()
    params.sim <- list()
    for (i in 1:nb.com){
      stats[[i]] <- lapply(models, function(x) x$sum.stats[[i]])
      params.sim[[i]] <- lapply(models, function(x) x$param)
    }
    stats <- lapply(stats, function(x) t(data.frame(x)))
    stats<- lapply(stats, function(x) {
      rownames(x) <- NULL
      return(x)
    })
    
    params.sim <- lapply(params.sim, function(x) t(data.frame(x)))
    params.sim <- lapply(params.sim, function(x) {
      rownames(x) <- NULL
      return(x)
    })
    params.sim <- lapply(params.sim, function(x, prior) {
      colnames(x) <- names(prior)
    return(x)
    }, prior = prior)
  }
  else {
    stats <- t(data.frame(lapply(models, function(x) x$sum.stats)))
    rownames(stats) <- NULL
    params.sim <- t(data.frame(lapply(models, function(x) x$param)))
    rownames(params.sim) <- NULL
    colnames(params.sim) <- names(prior)
  }
  
  return(list(stats = stats, params.sim = params.sim))
}
