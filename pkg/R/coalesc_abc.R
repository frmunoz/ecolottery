coalesc_abc <- function(comm.obs, pool, multi = F, traits=NULL, f.sumstats,
                        filt.abc, params, nb.samp = 10^6, parallel = T,
                        tol = 1*10^-4, pkg = NULL, method="neuralnet")
{
  
  if(is.null(pool)) {
    stop("You must provide regional pool composition")
  }
  
  if (!requireNamespace("abc", quietly = TRUE)) {
    stop("coalesc_abc requires package abc to be installed")
  }
  
  # Other required packages
  for (i in 1:length(pkg)) {
    
    if (!requireNamespace(pkg[i], quietly = TRUE)) {
      stop(paste("Package ", pkg[i], " is not available", sep = ""))
    }
    
  }
  
  # Community size
  if (!multi) {
    J <- nrow(comm.obs)
  } else if (var(tapply(comm.obs[, 2], comm.obs[, 1], length)) == 0) {
    J <- mean(tapply(comm.obs[, 2], comm.obs[, 1], length))
  } else {
    stop("multi option available with equal community sizes")
  }
  
  # Mean species traits
  if (ncol(pool) >= 3) {
    traits <- data.frame(apply(data.frame(pool[, -(1:2)]), 2,
                               function(x) tapply(x, pool[, 2], mean)))
  } else if (is.null(traits) & ncol(pool) < 3) {
    warning("Trait information is not provided")
  }
  
  # Compute summary statistics on multiple communities
  if (multi) {
    stats.obs <- c()
    for (i in unique(comm.obs[,1])) {
      if (is.null(traits)) {
        stats.obs <- rbind(stats.obs, f.sumstats(comm.obs[comm.obs[,1] == i,3]))
      } else {
        stats.obs <- rbind(stats.obs, f.sumstats(comm.obs[comm.obs[,1] == i,3],
                                                 traits[comm.obs[comm.obs[,1] == i,3],]))
      }
    }
  } else if (is.null(traits)) {
    stats.obs <- f.sumstats(comm.obs)
  } else {
    stats.obs <- f.sumstats(comm.obs[,2],traits[comm.obs[,2],])
  }
  
  # Community simulation
  sim <- do.simul(J, pool, traits, f.sumstats, filt.abc, params, nb.samp,
                  parallel, tol, pkg, method)
   
  # ABC estimation
  colnames(sim$stats) <- names(stats.obs)
  sim$stats.scaled <- sapply(1:ncol(sim$stats),
                             function(y) sim$stats[,y]/max(abs(sim$stats[,y]),
                                                           na.rm = T))
  colnames(sim$stats.scaled) <- colnames(sim$stats)
  meta.abc <- c()
  
  if (multi) {
    for (i in unique(comm.obs[,1])) {
    
      stats.obs.scaled <- sapply(1:length(stats.obs[i,]),
                                 function(x) stats.obs[i,x]/max(abs(sim$stats[,x]),
                                                                na.rm = T))
      meta.abc[[i]] <- abc::abc(target = stats.obs.scaled,
                                param = sim$params.sim,
                                sumstat = sim$stats.scaled,
                                tol = max(1.1/nrow(sim$stats.scaled), tol),
                                method = method)
    }
  } else {
    stats.obs.scaled <- sapply(1:length(stats.obs),
                               function(x) stats.obs[x]/max(abs(sim$stats[,x]),
                                                            na.rm = T))
    meta.abc <- abc::abc(target = stats.obs.scaled,
                         param = sim$params.sim,
                         sumstat = sim$stats.scaled,
                         tol = max(1.1/nrow(sim$stats.scaled), tol),
                         method = method)
  }
  
  return(list(par = sim$params.sim, ss = sim$stats.scaled, abc = meta.abc))
}

do.simul <- function(J, pool, traits = NULL, f.sumstats, filt.abc, params,
                     nb.samp = 10^6, parallel = T, tol = 1*10^-4, pkg = NULL,
                     method = "neuralnet") {
  if (!requireNamespace("parallel", quietly = TRUE) & parallel)  
  {
    warning("parallel = T requires package parallel to be installed; change to parallel = F")
    parallel <- F
  }
  
  if (parallel) {
    # Start up a parallel cluster
    parallelCluster <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
  }
  
  # Uniform prior distributions of parameters
  prior <- c()
  for (i in 1:nrow(params)) {
    prior[[i]] = runif(nb.samp, min = params[i, 1], max = params[i, 2])
  }              
  
  names(prior) <- rownames(params)
  prior[[length(prior) + 1]] <- runif(nb.samp, min = 0, max = 1)
  names(prior)[nrow(params) + 1] <- "m"
  
  # Function to perform simulations
  mkWorker <- function(traits, prior, J, pool, filt.abc, f.sumstats, pkg) {
    force(J)
    force(pool)
    force(traits)
    force(filt.abc)
    force(f.sumstats)
    force(prior)
    force(pkg)
    summCalc <- function(j,traits,prior,J,pool,filt.abc,f.sumstats) {
      params.samp <- unlist(lapply(prior,function(x) x[j]))
      stats.samp <- NA
      
      try({
        comm.samp <- coalesc(J, m = params.samp[length(params.samp)],
                             filt = function(x) filt.abc(x, params.samp),
                             pool = pool, traits = traits)
        stats.samp <- f.sumstats(comm.samp$com[, 2], comm.samp$com[, 3])
      })
      
      return(list(sum.stats = stats.samp, param = params.samp))
    }
    
    worker <- function(j) {
      require(lottery)
      # Other required packages
      for (i in 1:length(pkg)) {
        if (!requireNamespace(pkg[i], quietly = TRUE)) {
          stop(paste("Package ", pkg[i], " is not available", sep = ""))
          }
        }
      summCalc(j, traits, prior, J, pool, filt.abc, f.sumstats)
    }
    return(worker)
  }
  
  # Calculation of summary statistics over the whole range of parameters
  if (parallel) {
    models <- parallel::parLapply(parallelCluster, 1:nb.samp,
                                  mkWorker(traits, prior, J, pool, filt.abc,
                                           f.sumstats, pkg))
  } else {
    models <- lapply(1:nb.samp, mkWorker(traits, prior, J, pool, filt.abc,
                                         f.sumstats, pkg))
  }
  
  if (parallel) {
    # Close parallelCluster
    if (!is.null(parallelCluster)) {
      parallel::stopCluster(parallelCluster)
      parallelCluster <- c()
    }
  }
  
  stats <- t(data.frame(lapply(models, function(x) x$sum.stats)))
  rownames(stats) <- NULL
  params.sim <- t(data.frame(lapply(models, function(x) x$param)))
  rownames(params.sim) <- NULL
  colnames(params.sim) <- c(row.names(params), "m")
  
  return(list(stats = stats, params.sim = params.sim))
}
