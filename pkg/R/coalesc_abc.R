coalesc_abc <- function(comm.obs, pool = NULL, multi = "none", traits = NULL,
                        f.sumstats, filt.abc = NULL, params, theta.max = NULL,
                        nb.samp = 10^6, parallel = TRUE, tol = 1*10^-4, 
                        pkg = NULL, method = "neuralnet")
{
  
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
  if (!(multi %in% c("none", "tab", "seqcom"))){
    stop("multi parameter must be set at none, tab or seqcom.")
  }
  
  if (multi == "none") {
    J <- nrow(comm.obs)
    nb.com <- 1
  } else {
    if (multi == "tab") { 
      if (ncol(comm.obs) < 3){
        stop("When multi is set at tab, comm.obs must be a data.frame with at least three columns,
             including number of communities, number of individuals, number of species and optionally
             trait values.")
      }
      
      if (var(table(comm.obs[, 1])) == 0) {
        J <- mean(table(comm.obs[, 1]))
      } else {
        stop("multi option available with equal community sizes")
      }
      # nb.com is equal to number communities
      nb.com <- length(table(comm.obs[, 1]))
      # Deletion of first column with community number, so the third one contains trait values
      comm.obs <- comm.obs[, 2:ncol(comm.obs)]
    } else if (multi == "seqcom") stop("seqcom option not implemented yet")
  }
  
  # Trait values can be provided with community composition
  # Mean trait values of pool are stored in traits in absence of
  # trait information in local community
  if (!is.null(pool) & ncol(comm.obs) < 3 & ncol(pool) >= 3) {
      traits <- data.frame(apply(data.frame(pool[, -(1:2)]), 2,
                                 function(x) tapply(x, pool[, 2], mean)))
  }
  if (is.null(traits) & ncol(comm.obs) < 3) {
    warning("Trait information is not provided")
  }
  
  if (is.null(traits)) {
    stats.obs <- f.sumstats(comm.obs)
  } else {
    stats.obs <- f.sumstats(comm.obs, traits)
  }
  
  # Community simulation
  sim <- do.simul(J, pool, multi, nb.com, traits, f.sumstats, filt.abc, params,
                  theta.max, nb.samp, parallel, tol, pkg, method)
   
  stats.mean <- apply(sim$stats, 2, function(x) mean(x, na.rm = T))
  stats.sd <- apply(sim$stats, 2, function(x) sd(x, na.rm = T))
  sim$stats.scaled <- t(apply(sim$stats, 1,
                              function(x) (x - stats.mean)/stats.sd))
  
  colnames(sim$stats.scaled) <- colnames(sim$stats)
 
  stats.obs.scaled <- (stats.obs - stats.mean)/stats.sd
  
  # ABC estimation
  sel <- which(rowSums(is.na(sim$stats.scaled)) == 0)
  res.abc <- abc::abc(target = stats.obs.scaled,
                      param = sim$params.sim[sel,],
                      sumstat = sim$stats.scaled[sel,],
                      tol = tol,
                      method = method)
    
  return(list(par = sim$params.sim, ss = sim$stats.scaled, abc = res.abc))
}

do.simul <- function(J, pool = NULL, multi = "none", nb.com = NULL,
                     traits = NULL, f.sumstats = NULL, filt.abc = NULL, params,
                     theta.max = NULL, nb.samp = 10^6, parallel = TRUE,
                     tol = 1*10^-4, pkg = NULL, method = "neuralnet") {
  
  if (!requireNamespace("parallel", quietly = TRUE) & parallel) {
    warning(paste0("parallel = TRUE requires package 'parallel' to be ",
                   " installed\nchanged to parallel = FALSE"))
    parallel <- FALSE
  }
  
  if (parallel) {
    # Start up a parallel cluster
    parCluster <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
  }
  
  # Uniform prior distributions of parameters
  prior <- c()
  for (i in 1:nrow(params)) {
    prior[[i]] <- runif(nb.samp, min = params[i, 1], max = params[i, 2])
  }              
  names(prior) <- rownames(params)
  
  prior[[length(prior) + 1]] <- runif(nb.samp, min = 0, max = 1)
  names(prior)[nrow(params) + 1] <- "m"
  
  if (is.null(pool)) {
    prior[[length(prior) + 1]] <- runif(nb.samp, min = 0, max = theta.max)
    names(prior)[length(prior)] <- "theta"
  }
  
  # Function to perform simulations
  mkWorker <- function(traits, nb.com, multi, prior, J, pool, filt.abc,
                       f.sumstats, pkg) {
    force(J)
    force(pool)
    force(traits)
    force(filt.abc)
    force(f.sumstats)
    force(prior)
    force(nb.com)
    force(pkg)
    
    summCalc <- function(j, multi, traits, nb.com, prior, J, pool, filt.abc,
                         f.sumstats) {
      
      params.samp <- unlist(lapply(prior,function(x) x[j]))
      stats.samp <- NA
      params.samp.all <- params.samp
      names(params.samp.all) <- names(prior)
        
      if (is.null(pool)) {
        pool <- coalesc(J*100, theta = params.samp[length(params.samp)])$pool
        params.samp <- params.samp[-length(params.samp)]
      }
      
      if (!is.null(filt.abc)) {
        filt <- function(x) filt.abc(x, params.samp[-length(params.samp)])
      } else {
        filt <- NULL
      }
      
      if (nb.com > 1) {
        
        pool.sp <- unique(pool$sp)
        meta.samp <- array(0, c(nb.com, length(pool.sp)))
        colnames(meta.samp) <- pool.sp
                
        for (i in 1:nb.com) {
          
          try({
            comm.samp <- coalesc(J, m = params.samp[length(params.samp)],
                                 filt = filt,
                                 pool = pool, traits = traits)
            tab <- table(comm.samp$com[,2])
            meta.samp[i,names(tab)] <- tab
          })
        }
        
        if (is.null(traits)) {
          stats.samp <- f.sumstats(meta.samp)
        } else {
          stats.samp <- f.sumstats(meta.samp, traits)
        }
      
      } else {
        comm.samp <- coalesc(J, m = params.samp[length(params.samp)],
                             filt = filt,
                             pool = pool, traits = traits)
        if (is.null(traits)) {
          stats.samp <- f.sumstats(comm.samp$com)
        } else {
          stats.samp <- f.sumstats(comm.samp$com, traits)
        }
      }
          
      return(list(sum.stats = stats.samp, param = params.samp.all))
    }
    
    worker <- function(j) {
      require(lottery)
      # Other required packages
      if (!is.null(pkg)) for (i in 1:length(pkg)) {
        if (!requireNamespace(pkg[i], quietly = TRUE)) {
          stop(paste("Package ", pkg[i], " is not available", sep = ""))
        }
      }
      summCalc(j, multi, traits, nb.com, prior, J, pool, filt.abc, f.sumstats)
    }
    return(worker)
  }
  
  # Calculation of summary statistics over the whole range of parameters
  if (parallel) {
    models <- parallel::parLapply(parCluster, 1:nb.samp,
                                  mkWorker(traits, nb.com, multi, prior, J,
                                           pool, filt.abc, f.sumstats, pkg))
  } else {
    models <- lapply(1:nb.samp, mkWorker(traits, nb.com, multi, prior, J, pool,
                                         filt.abc, f.sumstats, pkg))
  }
  
  if (parallel) {
    # Close parCluster
    if (!is.null(parCluster)) {
      parallel::stopCluster(parCluster)
      parCluster <- c()
    }
  }
  
  stats <- t(data.frame(lapply(models, function(x) x$sum.stats)))
  rownames(stats) <- NULL
  params.sim <- t(data.frame(lapply(models, function(x) x$param)))
  rownames(params.sim) <- NULL
  colnames(params.sim) <- names(prior)
  
  return(list(stats = stats, params.sim = params.sim))
}
