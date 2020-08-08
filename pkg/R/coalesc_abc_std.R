coalesc_abc_std <- function(comm.obs, pool = NULL, multi = "single", prop = F, traits = NULL,
                            f.sumstats, filt.abc = NULL, filt.vect = F, migr.abc = NULL, size.abc = NULL, 
                            add = F, var.add = NULL, params = NULL, par.filt = NULL, par.migr = NULL, 
                            par.size = NULL, constr = NULL, scale = F, dim.pca = NULL, svd = F, theta.max = NULL, 
                            nb.samp = 10^6, parallel = TRUE, nb.core = NULL, tol = NULL, pkg = NULL, 
                            method.abc = "rejection")
{
  if(is.null(par.filt) & is.null(par.migr)){
    warning("No value provided for par.filt and par.migr arguments. ",
            "Only m and theta will be estimated.", call. = FALSE)
  }
  
  # The use of pool.glob must be checked
  if (!is.data.frame(pool) & is.list(pool) & length(pool)>1) {
    pool.glob <- Reduce(rbind, pool)
  } else {
    pool.glob <- pool
  }
  
  if (!is.null(pool.glob) & !is.null(traits) &
      sum(!pool.glob[,2] %in% rownames(traits)!=0)) {
    warning("The names of some species of the pool(s) are not present in ",
            "the rownames of traits", call. = FALSE)
  }
  
  # Community size
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
        stop("Species abundance must be integer values. ",
             "Consider using prop = TRUE for proportion data", call. = FALSE)
      nb.com <- nrow(comm.obs)
    } else if (multi == "seqcom") {
      # comm.obs is a list of communities with individuals on rows in each
      # community
      J <- unlist(lapply(comm.obs, nrow))
      if(any(round(J) != J)) 
        stop("Species abundance must be integer values. ",
             "Consider using prop = TRUE for proportion data", call. = FALSE)
      nb.com <- length(comm.obs)
    }
  }
  
  if (!is.null(pool.glob) & is.null(traits))
    # Takes average trait values for species defined in pool
    if(ncol(pool.glob) > 3) {
      traits <- data.frame(apply(pool.glob[,-(1:2)], 2,
                                 function(x) {
                                   tapply(x, pool.glob[, 2],
                                          function(y)
                                            mean(y, na.rm = TRUE)
                                   )}), row.names = unique(pool.glob[, 2]))
    } else if(ncol(pool.glob) == 3) {
      traits <- data.frame(tapply(pool.glob[,-(1:2)], pool.glob[, 2],
                                  function(y)
                                    mean(y, na.rm = TRUE)), row.names = unique(pool.glob[, 2]))
    } 
  
  if (length(formals(f.sumstats)) == 1) {
    stats.obs <- f.sumstats(comm.obs)
  } else if (length(formals(f.sumstats)) == 2) {
    stats.obs <- f.sumstats(comm.obs, traits)
  } else {
    stats.obs <- f.sumstats(comm.obs, traits, var.add)
  }
  
  # Community simulation
  sim <- do.simul.coalesc(J, pool, multi, prop, nb.com, traits, f.sumstats, filt.abc, filt.vect, migr.abc,
                          size.abc, add, var.add, params, par.filt, par.migr, par.size, constr,
                          dim.pca, svd, theta.max, nb.samp, parallel, nb.core, tol, pkg, method.abc)
  
  if(sum(sim$sel.ss)!=length(stats.obs))
  {
    # Remove summary statistics that failed in simulation
    stats.obs <- stats.obs[sim$sel.ss]
    
    warning("Some summary statistics yielded many NA values and have been ",
            "withdrawn. Please consider redifining these statistics or ",
            "changing the prior distributions.", call. = FALSE)
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
        if(nrow(sim$stats.scaled)!=nrow(sim$params.sim)) {
        stop("stats.scaled and params.sim must have the same number of rows",
             call. = FALSE)
        }
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
                 method = method.abc),
        error = function(x) {
          warning("ABC computation failed with the requested method.",
                  call. = FALSE)
        })
    } else {
      res.abc <- tryCatch(
        abc::abc(target = stats.obs,
                 param = sim$params.sim,
                 sumstat = sim$stats,
                 tol = tol,
                 method = method.abc),
        error = function(x) {
          warning("ABC computation failed with the requested method.",
                  call. = FALSE)
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
                             traits = NULL, f.sumstats = NULL, filt.abc = NULL, filt.vect = F, migr.abc = NULL, 
                             size.abc = NULL, add = F, var.add = NULL, params = NULL, par.filt = NULL, 
                             par.migr = NULL, par.size = NULL, constr = NULL, dim.pca = NULL, svd = F, 
                             theta.max = NULL, nb.samp = 10^6, parallel = TRUE, nb.core = NULL, tol = NULL, 
                             pkg = NULL, method.abc = "rejection") 
{
  
  if (!requireNamespace("parallel", quietly = TRUE) & parallel) {
    warning("parallel = TRUE requires package 'parallel' to be installed\n",
            "changed to parallel = FALSE", call. = FALSE)
    parallel <- FALSE
  }
  
  # This part is problematic if the user wants that the second argument of f.sumstats
  # incorporates intraspecific trait variation
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
      stop("Missing trait information", call. = FALSE)
    }
  }
  
  if(parallel) {
    # Start up a parallel cluster
    if(is.null(nb.core)) nb.core <- parallel::detectCores() - 1
    parCluster <- parallel::makeCluster(max(1, nb.core))
  }
  
  # Generate a set of parameter values drawn from prior distributions
  prior <- generate_prior(pool, prop, constr, params, par.filt, par.migr, par.size, 
                          theta.max, nb.samp)
  
  # Function to perform simulations
  mkWorker <- function(traits, nb.com, multi, prop, prior, J, pool, filt.abc, filt.vect, 
                       migr.abc, size.abc, add, var.add, f.sumstats, pkg) {
    force(traits)
    force(nb.com)
    force(multi)
    force(prop)
    force(prior)
    force(J)
    force(pool)
    force(filt.abc)
    force(filt.vect)
    force(migr.abc)
    force(size.abc)
    force(add)
    force(var.add)
    force(f.sumstats)
    force(pkg)
    
    summCalc <- function(j, multi, traits, nb.com, prior, J, prop, pool, filt.abc, filt.vect, migr.abc, add, var.add,
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
        if(multi == "seqcom") {
          J <- lapply(1:nb.com,function(x) J)
        }
      }
      
      if (is.null(pool)) {
        if (multi == "seqcom"){
          pool <- coalesc(mean(unlist(J))*100,
                          theta = params.samp[length(params.samp)], checks = F)$com
        } else {
          pool <- coalesc(mean(J)*100,
                          theta = params.samp[length(params.samp)], checks = F)$com
        }
        params.samp <- params.samp[-length(params.samp)]
      }
      
     if(!is.null(filt.abc))  filt <- ifelse(!add, 
                            function(x) filt.abc(x, params.samp[1:nrow(par.filt)]),
                            function(x, var.add) filt.abc(x, params.samp[1:nrow(par.filt)], var.add))
     else filt <- NULL
      migr <- ifelse(!is.null(migr.abc),
                     ifelse(!add, 
                            function() migr.abc(params.samp[(nrow(par.filt)+1):(nrow(par.filt)+nrow(par.migr))]),
                            function(var.add) migr.abc(params.samp[(nrow(par.filt)+1):(nrow(par.filt)+nrow(par.migr))], var.add)),
                     ifelse(!add,
                            function() params.samp["m"],
                            function(var.add) params.samp["m"]))
      size <- ifelse(!is.null(size.abc),
                     ifelse(!add, 
                            function() size.abc(params.samp[(nrow(par.filt)+1):(nrow(par.filt)+nrow(par.migr))]),
                            function(var.add) size.abc(params.samp[(nrow(par.filt)+1):(nrow(par.filt)+nrow(par.migr))], var.add)),
                     ifelse(!add,
                            function() function() params.samp["J"],
                            function(var.add) params.samp["J"]))
      
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
                                             filt.vect = filt.vect,
                                             add = add,
                                             var.add = var.add[i,],
                                             pool = pool.loc, traits = traits,
                                             checks = F)
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
            if(!is.data.frame(pool) & is.list(pool) & length(pool)>1) {
              pool.loc <- pool[[i]]
            } else {
              pool.loc <- pool
            }
            seqcom.samp[[i]] <- ecolottery::coalesc(J[[i]],
                                                    m = m,
                                                    filt = filt,
                                                    filt.vect = filt.vect,
                                                    add = add,
                                                    var.add = var.add[i,],
                                                    pool = pool.loc, traits = traits, checks = F)$com
            if(prop) {
              # The result is a table with relative species proportions and trait values
              tr <- tapply(seqcom.samp[[i]][,3], seqcom.samp[[i]][,2], mean)
              seqcom.samp[[i]] <- data.frame(sp=names(tr), 
                                             cov=tapply(seqcom.samp[[i]][,3], seqcom.samp[[i]][,2], length)/J[[i]],
                                             tr=tr)
            }
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
                                         filt.vect = filt.vect,
                                         add = add,
                                         var.add = var.add,
                                         pool = pool, traits = traits,
                                         checks = F)
        if(prop) {
          # The result is a table with relative species proportions and trait values
          tr <- tapply(comm.samp$com[,3], comm.samp$com[,2], mean)
          comm.samp$com <- data.frame(sp=names(tr), 
                                         cov=tapply(comm.samp$com[,3], comm.samp$com[,2], length)/J[[i]],
                                         tr=tr)
        }
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
        stop("Package ecolottery is not available", call. = FALSE)
      }
      # Other required packages
      if (!is.null(pkg)) for (i in 1:length(pkg)) {
        if (!requireNamespace(pkg[i], quietly = TRUE)) {
          stop(paste("Package ", pkg[i], " is not available", sep = ""),
               call. = FALSE)
        }
      }
      summCalc(j, multi, traits, nb.com, prior, J, prop, pool, filt.abc, filt.vect, migr.abc, 
               add, var.add, f.sumstats)
    }
    return(worker)
  }
  
  # Calculation of summary statistics over the whole range of parameters
  if (parallel) {
    err.chk <- try(models <- parallel::parLapply(parCluster, 1:nb.samp,
                                                 mkWorker(traits, nb.com, multi, prop, prior, J,
                                                          pool, filt.abc, filt.vect, migr.abc, size.abc, add, var.add, f.sumstats, pkg)), T)
  } else {
    models <- lapply(1:nb.samp, mkWorker(traits, nb.com, multi, prop, prior, J, pool,
                                         filt.abc, filt.vect, migr.abc, size.abc, add, var.add, f.sumstats, pkg))
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