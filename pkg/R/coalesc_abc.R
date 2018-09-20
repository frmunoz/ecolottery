coalesc_abc <- function(comm.obs, pool = NULL, multi = "single", prop = F, traits = NULL,
                        f.sumstats, filt.abc = NULL, migr.abc = NULL, size.abc = NULL, add = F,
                        var.add = NULL, params = NULL, par.filt = NULL, par.migr = NULL, 
                        par.size = NULL, constr = NULL, scale = T, dim.pca = NULL, svd = F, theta.max = NULL, 
                        nb.samp = 10^6, parallel = TRUE, nb.core = NULL, tol = NULL, pkg = NULL, 
                        method = "rejection")
{
  if(!method%in%c("rejection", "loclinear", "neuralnet", "ridge"))
    stop("method.abc should be either rejection, loclinear, neuralnet or ridge")
  
  if(!is.function(f.sumstats)) {
    stop("You must provide a function to calculate summary statistics",
         "(f.sumstats)")
  }
  
  if(prop & multi!="tab") {
    stop("prop data can only be handled in tab format")
  }
  
  if(multi=="tab") if(any(rowSums(comm.obs)==0)) {
    stop("There should not be communities with 0 individuals")
  }
  if (length(formals(f.sumstats)) > 3) {
    stop("f.sumstats must be a function of up to three arguments")
  }
  
  if (is.null(tol)){
    warning("You must provide a tolerance value for ABC computation.\n",
            "The function will only provide simulations and will not perform ",
            "ABC analysis.")
  } 
  
  if(is.null(par.filt) & is.null(par.migr)){
    warning("No value provided for par.filt and par.migr arguments. Only m and theta will be ",
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
  
  if (!multi == "single" & class(pool) == "list") {
    if((multi == "tab" & length(pool) != nrow(comm.obs)) | (multi == "seqcom" & length(pool) != length(comm.obs)))
      stop("When several pools are given (list), there must be one for each community")
  }
  if (multi=="single" & class(pool)=="list" & length(pool)>1)
    stop("If multi=single, pool should not be a list of more than one element")
  
  if (!is.data.frame(pool) & is.list(pool) & length(pool)>1) {
    pool.glob <- Reduce(rbind, pool)
  } else {
    pool.glob <- pool
  }
  
  if (!is.null(pool.glob) & !is.null(traits) & sum(!pool.glob[,2]%in%rownames(traits)!=0)) {
    warning("The names of some species of the pool(s) are not present in the rownames of traits")
  }
  
  if (!requireNamespace("abc", quietly = TRUE)) {
    stop("coalesc_abc requires package abc to be installed")
  }
  
  if (!is.null(constr) & !requireNamespace("lazyeval", quietly = TRUE)) {
    stop("coalesc_abc requires package lazy_eval to be installed")
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
  if (!is.null(pool.glob)){
    if(ncol(pool.glob) >= 3) {
      #Using matrix instead of data.frame
      traits <- apply(data.frame(pool.glob[,-(1:2)]), 2,
                                 function(x) {
                                   tapply(x, pool.glob[, 2],
                                          function(y)
                                            mean(y, na.rm = TRUE)
                                   )})
    }
  }
  
  if (multi == "tab"){
    if(is.null(colnames(comm.obs)))
      stop("Missing species names in species-by-site table")
    if(any(!colnames(comm.obs)%in%rownames(traits)))
       stop("Mismatch of species names in pool and comm.obs")
    # Reorder species in comm.obs following the order in traits
    comm.obs <- comm.obs[,rownames(traits)[rownames(traits)%in%colnames(comm.obs)]]
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
  sim <- do.simul.coalesc(J, pool, multi, prop, nb.com, traits, f.sumstats, filt.abc, migr.abc,
                          size.abc, add, var.add, params, par.filt, par.migr, par.size, constr,
                          dim.pca, svd, theta.max, nb.samp, parallel, nb.core, tol, pkg, method)
  
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
  
  if (scale)
  {
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
  }
    
  if (is.null(tol)){
    res.abc <- NA
  } else {
    if(scale) {
      # ABC estimation
      res.abc <- tryCatch(
        abc::abc(target = stats.obs.scaled,
               param = sim$params.sim,
               sumstat = sim$stats.scaled,
               tol = tol,
               method = method),
        error = function(x) {
          warning("ABC computation failed with the requested method.")
        })
    } else {
      res.abc <- tryCatch(
        abc::abc(target = stats.obs,
                 param = sim$params.sim,
                 sumstat = sim$stats,
                 tol = tol,
                 method = method),
        error = function(x) {
          warning("ABC computation failed with the requested method.")
        })
    }
    
    if (is.character(res.abc)){
      res.abc <- NA
    }
  }
  
  if(is.null(dim.pca))
  {
    if(scale) {
      return(list(par = sim$params.sim, obs = stats.obs,
                obs.scaled = stats.obs.scaled, ss = sim$stats.scaled,
                ss.scale = data.frame(mean=stats.mean,sd=stats.sd),
                abc = res.abc, call = match.call()))
    } else {
      return(list(par = sim$params.sim, obs = stats.obs,
                  ss = sim$stats, abc = res.abc, call = match.call()))
    }
  } else
  {
    if(scale) {
      return(list(par = sim$params.sim, obs = stats.obs,
                obs.scaled = stats.obs.scaled, ss = sim$stats.scaled,
                abc = res.abc, call = match.call()))
    } else {
      return(list(par = sim$params.sim, obs = stats.obs,
                  ss = sim$stats, abc = res.abc, call = match.call()))
    }
  }
}

do.simul.coalesc <- function(J, pool = NULL, multi = "single", prop = F, nb.com = NULL,
                     traits = NULL, f.sumstats = NULL, filt.abc = NULL, migr.abc = NULL, 
                     size.abc = NULL, add = F, var.add = NULL, params = NULL, par.filt = NULL, 
                     par.migr = NULL, par.size = NULL, constr = NULL, dim.pca = NULL, svd = F, 
                     theta.max = NULL, nb.samp = 10^6, parallel = TRUE, nb.core = NULL, tol = NULL, 
                     pkg = NULL, method = "rejection") {
  
  if (!requireNamespace("parallel", quietly = TRUE) & parallel) {
    warning("parallel = TRUE requires package 'parallel' to be installed\n",
            "changed to parallel = FALSE")
    parallel <- FALSE
  }
  
  if (!is.null(params) & is.null(par.filt))
  {
    warning("The use of params will be deprecated, please use arguments par.filt 
            and par.migr instead.")
    par.filt <- params
    par.migr <- t(data.frame(c(0, 1)))
    rownames(par.migr) <- "m"
  }
  if (is.null(params) & is.null(par.filt) & is.null(par.migr) & (!is.null(filt.abc) | !is.null(migr.abc)))
  {
    stop("You must provide range limit values of the parameters.")
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
  
  # Defining constraints on parameter values
  constr.test <- function(par.names, par.val) 
  {
    res <- T
    if(!is.null(constr)) {
      for(j in 1:length(par.names))
        assign(par.names[j], par.val[j], envir=environment())
      for(j in 1:length(constr)) 
        res <- res & lazyeval::lazy_eval(constr[j], environment())
    }
    return(res)
  }
  
  # Uniform prior distributions of parameters
  prior <- list()
  dim.prior <- ifelse(!is.null(nrow(par.filt)),nrow(par.filt),0) + 
    max(ifelse(!is.null(par.migr),nrow(par.migr),0),1) + 
    prop*max(ifelse(!is.null(par.size),nrow(par.size),0),1) + as.numeric(is.null(pool))
  length(prior) <- dim.prior
  stop <- 0
  while(length(prior[[1]]) < nb.samp | stop < 100) {
    samp <- nb.samp - length(prior[[1]])
    i <- 0
    if(!is.null(par.filt)) {
      for (i in 1:nrow(par.filt)) {
        prior[[i]] <- c(prior[[i]], runif(samp, min = par.filt[i, 1], max = par.filt[i, 2]))
      }
    }
    if(!is.null(par.migr)) {
      for (i in 1:nrow(par.migr)) {
        prior[[i+nrow(par.filt)]] <- c(prior[[i+nrow(par.filt)]], runif(samp, min = par.migr[i, 1], 
                                                                        max = par.migr[i, 2]))
      }
    } else {
      prior[[i+1]] <- c(prior[[i+1]], runif(samp, min = 0, max = 1))
    }
    if(!is.null(par.migr)) {
      names(prior) <- c(rownames(par.filt), rownames(par.migr))
    } else {
      names(prior) <- c(rownames(par.filt), "m")
    }
    
    # Note - Defining a lower bound at 0 for m and theta can entail issues when species richness is 1
    # in simulated community; we should allow the user to define the prior for m and theta in
    # the future
    
    l <- length(prior)
    if(prop) 
    {
      if(!is.null(par.size))
      {
        for (i in 1:nrow(par.size)) {
          prior[[i+l]] <- c(prior[[i+l]], runif(prop, min = par.size[i, 1], max = par.size[i, 2]))
        }
        names(prior)[l:(l+nrow(par.size))] <- rownames(par.size)
      } else {
        prior[[length(prior)+1]] <- c(prior[[length(prior)+1]], runif(prop, min = 100, max = 1000))
        warning("The prior of community size is uniform between 100 and 1000")
        names(prior)[length(prior)] <- "J"
      }
    }
    
    if (is.null(pool)) {
      prior[[length(prior) + 1]] <- c(prior[[length(prior) + 1]], runif(prop, min = 0, max = theta.max))
      names(prior)[length(prior)] <- "theta"
    }
    
    constr.sel <- sapply(1:nb.samp, function(x) constr.test(names(prior), unlist(lapply(prior,function(y) y[x]))))
    prior <- lapply(prior, function(x) x[constr.sel])
    
    stop <- stop + 1
  }
  
  if(stop==100) warning("Fail to meet constraints for nb.samp sets of parameter values")
  
  # Function to perform simulations
  mkWorker <- function(traits, nb.com, multi, prop, prior, J, pool, filt.abc, migr.abc, size.abc,
                       add, var.add, f.sumstats, pkg) {
    force(traits)
    force(nb.com)
    force(multi)
    force(prop)
    force(prior)
    force(J)
    force(pool)
    force(filt.abc)
    force(migr.abc)
    force(size.abc)
    force(add)
    force(var.add)
    force(f.sumstats)
    force(pkg)
    
    summCalc <- function(j, multi, traits, nb.com, prior, J, prop, pool, filt.abc, migr.abc, add, var.add,
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
        
        ## Need to remove reference to par.filt here
        if (!is.null(filt.abc)) {
          if(!add) 
            filt <- function(x) filt.abc(x, params.samp[1:nrow(par.filt)])
          else filt <- function(x, var.add) filt.abc(x, params.samp[1:nrow(par.filt)], var.add)
        } else {
          filt <- NULL
        }
        if (!is.null(migr.abc)) {
          if(!add) 
            migr <- function() migr.abc(params.samp[(nrow(par.filt)+1):(nrow(par.filt)+nrow(par.migr))])
          else migr <- function(var.add) migr.abc(params.samp[(nrow(par.filt)+1):(nrow(par.filt)+nrow(par.migr))], var.add)
        } else {
          if(!add) 
            migr <- function() params.samp["m"]
          else migr <- function(var.add) params.samp["m"]
        }
        if (!is.null(size.abc)) {
          if(!add) 
            size <- function() size.abc(params.samp[(nrow(par.filt)+1):(nrow(par.filt)+nrow(par.migr))])
          else migr <- function(var.add) size.abc(params.samp[(nrow(par.filt)+1):(nrow(par.filt)+nrow(par.migr))], var.add)
        } else {
          if(!add) 
            size <- function() params.samp["J"]
          else size <- function(var.add) params.samp["J"]
        }
        
        if(multi == "tab") {
          # It is too slow
          if(!is.data.frame(pool) & is.list(pool) & length(pool) > 1) {
            pool.sp <- unique(Reduce(rbind , pool)[,2])
          } else {
            pool.sp <- unique(pool[,2])
          }
          meta.samp <- array(0, c(nb.com, length(pool.sp)))
          colnames(meta.samp) <- pool.sp
          
          for (i in 1:nb.com) {
            try({
              J.loc <- ifelse(length(J)>1, J[i], J)
              m <- unlist(ifelse(add, migr(var.add[i,]), migr()))
              if(!is.data.frame(pool) & is.list(pool) & length(pool)>1) {
                pool.loc <- pool[[i]]
              } else {
                pool.loc <- pool
              }
              comm.samp <- ecolottery::coalesc(J.loc, m = m,
                                   filt = filt,
                                   add = add,
                                   var.add = var.add[i,],
                                  pool = pool.loc, traits = traits)
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
              m <- ifelse(add, migr(var.add[i,]), migr())[[1]]
              if(is.list(pool) & length(pool)>1) {
                pool.loc <- pool[[i]]
              } else {
                pool.loc <- pool
              }
              seqcom.samp[[i]] <- ecolottery::coalesc(J[[i]],
                                          m = m,
                                          filt = filt,
                                          add = add,
                                          var.add = var.add[i,],
                                          pool = pool.loc, traits = traits)$com
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
          m <- ifelse(add, migr(var.add[i,]), migr())[[1]]
          comm.samp <- ecolottery::coalesc(J,
                             m = m,
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
      summCalc(j, multi, traits, nb.com, prior, J, prop, pool, filt.abc, migr.abc, 
               add, var.add, f.sumstats)
    }
    return(worker)
  }
  
  # Calculation of summary statistics over the whole range of parameters
  if (parallel) {
    err.chk <- try(models <- parallel::parLapply(parCluster, 1:nb.samp,
                                  mkWorker(traits, nb.com, multi, prop, prior, J,
                                           pool, filt.abc, migr.abc, size.abc, add, var.add, f.sumstats, pkg)), T)
  } else {
    models <- lapply(1:nb.samp, mkWorker(traits, nb.com, multi, prop, prior, J, pool,
                                         filt.abc, migr.abc, size.abc, add, var.add, f.sumstats, pkg))
  }
  
  if (parallel) {
    # Close parCluster
    if (!is.null(parCluster)) {
      parallel::stopCluster(parCluster)
      parCluster <- c()
    }
    if(class(err.chk)=="try-error")
    {
      stop(err.chk[1])
    }
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
  stats.sel <- stats[sel.row, sel.ss, drop = FALSE]
  params.sim.sel  <- params.sim[sel.row, , drop = FALSE]
  
  return(list(stats = stats.sel, params.sim = params.sim.sel, sel.ss=sel.ss))
}
