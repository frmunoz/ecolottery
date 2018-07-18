coalesc_abc <- function(comm.obs, pool = NULL, multi = "single", prop = F, traits = NULL,
                        f.sumstats, filt.abc = NULL, add = F, var.add = NULL,
                        params = NULL, dim.pca = NULL, svd = F, theta.max = NULL, nb.samp = 10^6, 
                        parallel = TRUE, nb.core = NULL, tol = NULL, pkg = NULL, method = "rejection")
{
  
  if(!method%in%c("rejection", "loclinear", "neuralnet", "ridge"))
    stop("method.abc should be either rejection, loclinear, neuralnet or ridge")
  
  if (!is.function(f.sumstats)) {
    stop("You must provide a function to calculate summary statistics",
         "(f.sumstats)")
  }
  
  if(prop & multi!="tab") {
    stop("prop data can only be handled in tab format")
  }
  
  if(multi=="tab" & any(rowSums(comm.obs)==0)) {
    stop("There should not be communities with 0 individuals")
  }
  if (length(formals(f.sumstats)) > 3) {
    stop("f.sumstats must be a function of up to three arguments")
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
  
  if (!is.null(pool) & !is.null(traits) & sum(!pool[,2]%in%rownames(traits)!=0)) {
    warning("The names of some species of the pool are not present in the rownames of traits")
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
      #if(prop & any(J != 1)) 
      #  stop("Relative species abundances must sum to 1")
      if(!prop & any(round(J) != J)) 
        stop("Species abundance must be integer values. Consider using prop = T for proportion data")
      nb.com <- nrow(comm.obs)
    } else if (multi == "seqcom") {
      # comm.obs is a list of communities with individuals on rows in each
      # community
      J <- unlist(lapply(comm.obs, nrow))
      if(any(round(J) != J)) 
        stop("Species abundance must be integer values. Consider using prop = T for proportion data")
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
  
  if (length(formals(f.sumstats)) == 1) {
    stats.obs <- f.sumstats(comm.obs)
  } else if (length(formals(f.sumstats)) == 2) {
    stats.obs <- f.sumstats(comm.obs, traits)
  } else {
    stats.obs <- f.sumstats(comm.obs, traits, var.add)
  }
  
  # Community simulation
  sim <- do.simul.coalesc(J, pool, multi, prop, nb.com, traits, f.sumstats, filt.abc, add, var.add, 
                  params, dim.pca, svd, theta.max, nb.samp, parallel, nb.core, tol, pkg, method)
  
  if(sum(sim$sel.ss)!=length(stats.obs))
  {
    # Remove summary statistics that failed in simulation
    stats.obs <- stats.obs[sim$sel.ss]
    
    warning("Some summary statistics yielded many NA values and have been 
            withdrawn. Please consider redifining these statistics or changing the prior distributions.")
  }
  
  # Scaling f.sumstats criterions for simulations
  if(is.null(dim.pca))
  {
    stats.mean <- apply(sim$stats, 2, function(x) mean(x, na.rm = TRUE))
    stats.sd <- apply(sim$stats, 2, function(x) sd(x, na.rm = TRUE))
    sim$stats.scaled <- t(apply(sim$stats, 1,
                              function(x) (x - stats.mean)/stats.sd))
  } else 
  {
    if(!svd)
    {
        stats.pca <- ade4::dudi.pca(rbind(sim$stats, stats.obs), scannf = F, nf = dim.pca) 
        # Use scores on PCA dimensions
        sim$stats.scaled <- sim$stats.pca$l1
    } else
    {
      bigtab <- rbind(sim$stats, stats.obs)
      bigtab <- scale(bigtab)
      stats.svd <- svd(bigtab)
      # Use scores derived from SVD
      sim$stats.scaled <- (stats.svd$u%*%diag(stats.svd$d))[,1:dim.pca]
    }
  }
  colnames(sim$stats.scaled) <- colnames(sim$stats)
  
  if (is.null(tol)){
    res.abc <- NA
    stats.obs.scaled <- (stats.obs - stats.mean)/stats.sd
  } else {
    if(is.null(dim.pca))
    {
      stats.obs.scaled <- (stats.obs - stats.mean)/stats.sd
    } else
    {
      stats.obs.scaled <- sim$stats.scaled[nrow(sim$stats.scaled),]
      sim$stats.scaled <- sim$stats.scaled[-nrow(sim$stats.scaled),]
      # For debug
      if(nrow(sim$stats.scaled)!=nrow(sim$params.sim)) stop("stats.scaled and params.sim must have
                                                              the same number of rows")
    }
  }
    
  if (is.null(tol)){
    res.abc <- NA
  } else {
    # ABC estimation
    res.abc <- tryCatch(
      abc::abc(target = stats.obs.scaled,
               param = sim$params.sim,
               sumstat = sim$stats.scaled,
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
  
  if(is.null(dim.pca))
  {
    return(list(par = sim$params.sim, obs = stats.obs,
                obs.scaled = stats.obs.scaled, ss = sim$stats.scaled,
                ss.scale = data.frame(mean=stats.mean,sd=stats.sd),
                abc = res.abc))
  } else
  {
    return(list(par = sim$params.sim, obs = stats.obs,
                obs.scaled = stats.obs.scaled, ss = sim$stats.scaled,
                abc = res.abc))
  }
}

do.simul.coalesc <- function(J, pool = NULL, multi = "single", prop = F, nb.com = NULL,
                     traits = NULL, f.sumstats = NULL, filt.abc = NULL, 
                     add = F, var.add = NULL, params, dim.pca = NULL, svd = F,
                     theta.max = NULL, nb.samp = 10^6, parallel = TRUE, nb.core = NULL, 
                     tol = NULL, pkg = NULL, method = "rejection") {
  
  if (!requireNamespace("parallel", quietly = TRUE) & parallel) {
    warning("parallel = TRUE requires package 'parallel' to be installed\n",
            "changed to parallel = FALSE")
    parallel <- FALSE
  }
  
  if(length(formals(f.sumstats))>1 & is.null(traits))
  {
    if(ncol(pool) >= 3) {
      traits <- data.frame(apply(data.frame(pool[,-(1:2)]), 2,
                                 function(x) {
                                   tapply(x, pool[, 2],
                                          function(y)
                                            mean(y, na.rm = TRUE)
                                   )}))
    } else {
      stop("Missing trait information")
    }
  }
  
  if(parallel) {
    # Start up a parallel cluster
    if(is.null(nb.core)) nb.core <- parallel::detectCores() - 1
    parCluster <- parallel::makeCluster(max(1, nb.core))
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
  if(!is.null(params)) {
    names(prior)[nrow(params) + 1] <- "m"
  } else
  {
    names(prior) <- "m"
  }
  
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
  mkWorker <- function(traits, nb.com, multi, prop, prior, J, pool, filt.abc, add, var.add,
                       f.sumstats, pkg) {
    force(traits)
    force(nb.com)
    force(multi)
    force(prop)
    force(prior)
    force(J)
    force(pool)
    force(filt.abc)
    force(add)
    force(var.add)
    force(f.sumstats)
    force(pkg)
    
    summCalc <- function(j, multi, traits, nb.com, prior, J, prop, pool, filt.abc, add, var.add,
                         f.sumstats) {
      
        params.samp <- unlist(lapply(prior,function(x) x[j]))
        stats.samp <- NA
        params.samp.all <- params.samp
        names(params.samp.all) <- names(prior)
        
        if(prop) 
        {
          if(!is.null(pool))
          {
            J <- round(params.samp[length(params.samp)])
            params.samp <- params.samp[-length(params.samp)]
          } else 
          {
            J <- round(params.samp[length(params.samp)-1])
            params.samp <- params.samp[-(length(params.samp)-1)]
          }
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
          if(!add) 
            filt <- function(x) filt.abc(x, params.samp[-length(params.samp)])
          else filt <- function(x, var.add) filt.abc(x, params.samp[-length(params.samp)], var.add)
        } else {
          filt <- NULL
        }
        
        if(multi == "tab") {
          pool.sp <- unique(pool[,2])
          meta.samp <- array(0, c(nb.com, length(pool.sp)))
          colnames(meta.samp) <- pool.sp
          
          for (i in 1:nb.com) {
            try({
              J.loc <- ifelse(length(J)>1, J[i], J)
              comm.samp <- ecolottery::coalesc(J.loc, m = params.samp[length(params.samp)],
                                   filt = filt,
                                   add = add,
                                   var.add = var.add[i,],
                                  pool = pool, traits = traits)
              tab <- table(comm.samp$com[,2])
              meta.samp[i,names(tab)] <- tab
              if(prop) meta.samp[i,] <- meta.samp[i,]/J.loc #sum(meta.samp[i,])
            })
          }
          
          if (length(formals(f.sumstats)) == 1) {
            stats.samp <- f.sumstats(meta.samp)
          } else if (length(formals(f.sumstats)) == 2) {
            stats.samp <- f.sumstats(meta.samp, traits)
          } else {
            stats.samp <- f.sumstats(meta.samp, traits, var.add)
          }
          
        } else if(multi == "seqcom") {
          seqcom.samp <- list()
          
          for (i in 1:nb.com) {
            try({
              seqcom.samp[[i]] <- ecolottery::coalesc(J[[i]],
                                          m = params.samp[length(params.samp)],
                                          filt = filt,
                                          add = add,
                                          var.add = var.add[i,],
                                          pool = pool, traits = traits)$com
            })
          }
          if (length(formals(f.sumstats)) == 1) {
            stats.samp <- f.sumstats(seqcom.samp)
          } else if (length(formals(f.sumstats)) == 2) {
            stats.samp <- f.sumstats(seqcom.samp, traits)
          } else {
            stats.samp <- f.sumstats(seqcom.samp, traits, var.add)
          }
        } else { # single community
        comm.samp <- ecolottery::coalesc(J,
                             m = params.samp[length(params.samp)],
                             filt = filt,
                             add = add,
                             var.add = var.add,
                             pool = pool, traits = traits)
        if(prop) comm.samp$com <- t(table(comm.samp$com[,2])/J)
        if (length(formals(f.sumstats)) == 1) {
          stats.samp <- f.sumstats(comm.samp$com)
        } else if (length(formals(f.sumstats)) == 2) {
          stats.samp <- f.sumstats(comm.samp$com, traits)
        } else {
          stats.samp <- f.sumstats(comm.samp$com, traits, var.add)
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
      summCalc(j, multi, traits, nb.com, prior, J, prop, pool, filt.abc, 
               add, var.add, f.sumstats)
    }
    return(worker)
  }
  
  # Calculation of summary statistics over the whole range of parameters
  if (parallel) {
    err.chk <- try(models <- parallel::parLapply(parCluster, 1:nb.samp,
                                  mkWorker(traits, nb.com, multi, prop, prior, J,
                                           pool, filt.abc, add, var.add, f.sumstats, pkg)), T)
  } else {
    models <- lapply(1:nb.samp, mkWorker(traits, nb.com, multi, prop, prior, J, pool,
                                         filt.abc, add, var.add, f.sumstats, pkg))
  }
  
  if (parallel) {
    # Close parCluster
    if (!is.null(parCluster)) {
      parallel::stopCluster(parCluster)
      parCluster <- c()
    }
  }
  if(class(err.chk)=="try-error")
  {
    stop(err.chk[1])
  }
  
  stats <- t(data.frame(lapply(models, function(x) x$sum.stats)))
  rownames(stats) <- NULL
  params.sim <- t(data.frame(lapply(models, function(x) x$param)))
  rownames(params.sim) <- NULL
  colnames(params.sim) <- names(prior)
  
  # Remove simulations for which some summary statistics are NA
  #sel <- which(rowSums(is.na(stats)) == 0)
  #params.sim <- params.sim[sel,]
  #stats <- stats[sel,]
  # Remove summary statistics with more than 50% NA
  sel.ss <- colSums(is.na(stats))<nrow(stats)/2
  # Remove rows with NA
  sel.row <- rowSums(is.na(data.frame(stats[, sel.ss])))==0
  stats.sel <- stats[sel.row, sel.ss]
  params.sim.sel  <- params.sim[sel.row, ]
  
  return(list(stats = stats.sel, params.sim = params.sim.sel, sel.ss=sel.ss))
}
