coalesc_abc <- function(comm.obs, pool = NULL, multi = "single", prop = FALSE,
                        traits = NULL, f.sumstats, filt.abc = NULL,
                        migr.abc = NULL, size.abc = NULL, add = FALSE,
                        var.add = NULL, params  = NULL, par.filt = NULL,
                        par.migr = NULL, par.size = NULL, scale = FALSE,
                        dim.pca = NULL, svd = FALSE, theta.max = NULL,
                        nb.samp = 10^6, parallel = FALSE, nb.core = NULL,
                        tol = NULL, type = "standard", method.seq = "Lenormand",
                        method.mcmc = "Marjoram_original",
                        method.abc = "rejection", alpha = 0.5, pkg = NULL) 
{
  
  if (is.null(pool)) {
    warning("No species pool provided: pool will be simulated and logseries ",
            "theta estimated", call. = FALSE)
    if (is.null(theta.max))
      theta.max <- 500
  }
  
  if (is.character(comm.obs) & is.null(ncol(comm.obs))) {
    comm.obs <- data.frame(id = 1:length(comm.obs),
                           sp = comm.obs,
                           trait = rep(NA, length(comm.obs)),
                           stringsAsFactors = FALSE)
    colnames(comm.obs) <- c("id", "sp", "trait")
  }
  
  if (multi == "tab"){
    if(is.null(colnames(comm.obs)))
      warning("Missing species names in species-by-site table", call. = FALSE)
    if(length(pool) == 1) if(any(!colnames(comm.obs) %in% unique(pool$sp))) 
      stop("Mismatch of species names in pool and comm.obs", call. = FALSE)
    else if(any(!colnames(comm.obs) %in% unique(Reduce(rbind, pool)$sp)))
      stop("Mismatch of species names in pools and comm.obs", call. = FALSE)
    if(!is.null(traits))
    {
      if(any(!colnames(comm.obs)%in%row.names(traits)))
        stop("Mismatch of species names in comm.obs and traits", call. = FALSE)
      # Reorder species in comm.obs following the order in traits
      comm.obs <- comm.obs[,rownames(traits)[rownames(traits) %in%
                                               colnames(comm.obs)]]
    }
  }
  
  # Defining parameter ranges
  if (!is.null(params) & is.null(par.filt))
  {
    warning("The use of params will be deprecated, please use arguments ",
            "par.filt and par.migr instead.", call. = FALSE)
    par.filt <- params
  }
  if (is.null(par.migr))
  {
    par.migr <- t(data.frame(c(0, 1)))
    rownames(par.migr) <- "m"
  }
  
  if(!is.null(nb.core)) if(nb.core == 1) parallel <- F
  
  if(is.null(type)) {
    warning('Standard ABC analysis (type = "standard")', call. = FALSE)
    type <- "standard"
  }
  
  # The options migr.abc = NULL, size.abc = NULL, add = F, 
  # var.add = NULL, par.filt = NULL, par.migr = NULL, 
  # par.size = NULL, constr = NULL, scale = T, dim.pca = NULL, svd = F, pkg = NULL
  # need to be included
  if (type == "standard")
    return(coalesc_abc_std(comm.obs, pool, multi, prop, traits, f.sumstats,
                           filt.abc, migr.abc, size.abc, add, var.add, params,
                           par.filt, par.migr, par.size, constr = NULL,
                           scale = TRUE, dim.pca, svd, theta.max, nb.samp,
                           parallel, nb.core, tol, pkg, method.abc)) 
  else {
    if(!is.null(migr.abc) | !is.null(size.abc))
    {
      migr.abc <- NULL
      size.abc <- NULL
      warning("migr.abc, size.abc not available with other sampling approach ",
              "than standard", call. = FALSE)
    }
    if(add)
    {
      add <- FALSE
      warning("add and var.add not available with other sampling approach ",
              "than standard", call. = FALSE)
    }
    if(!is.null(dim.pca) | svd)
    {
      dim.pca <- NULL
      svd <- FALSE
      warning("dim.pca and svd not available with other sampling approach ",
              "than standard", call. = FALSE)
    } 
    initial_checks(comm.obs = comm.obs, pool = pool, multi = multi, prop = prop,
                   traits = traits, f.sumstats = f.sumstats,
                   filt.abc = filt.abc, migr.abc = migr.abc,
                   size.abc = size.abc, params = params, par.filt = par.filt,
                   par.migr = par.migr,  par.size = par.size, scale = scale,
                   theta.max = theta.max, nb.samp = nb.samp,
                   parallel = parallel, nb.core = nb.core, tol = tol,
                   type = type, method.seq = method.seq,
                   method.mcmc = method.mcmc, method.abc = method.abc,
                   alpha = alpha, pkg = pkg)
  }  
  
  if (parallel){
    stop("parallel computation only implemented for type = standard - ",
         "ongoing work", call. = FALSE)
  }
  
  if (is.null(pool)) {
    warning("No species pool provided: pool will be simulated and logseries ",
            "theta estimated", call. = FALSE)
    if (is.null(theta.max))
      theta.max <- 500
    stop("Estimation of theta not implemented - Ongoing work", call. = FALSE)
  }
  
  if(prop & multi!="tab")
    stop("prop data can only be handled in tab format", call. = FALSE)
  
  if (!requireNamespace("EasyABC", quietly = TRUE)) 
    stop("coalesc_abc requires package EasyABC to be installed", call. = FALSE)
  #for (i in 1:length(pkg)) if (!requireNamespace(pkg[i], quietly = TRUE)) 
  #  stop(paste("Package ", pkg[i], " is not available", sep = ""))
  
  if(!type %in% c("seq","mcmc","annealing","standard"))
    stop("type should be either standard, seq, mcmc or annealing",
         call. = FALSE)
  
  if(type == "seq") if(!method.seq %in% c("Beaumont", "Drovandi", "Delmoral",
                                          "Lenormand", "Emulation"))
    stop("method.seq should be either Beaumont, Drovandi, Delmoral, Lenormand ",
         "or Emulation", call. = FALSE)
  
  if(type == "mcmc") if(!method.mcmc %in% c("Marjoram_original", "Marjoram",
                                            "Wegmann"))
    stop("method.mcmc should be either Marjoram_original, Marjoram or Wegmann",
         call. = FALSE)
  
  if(!is.null(nb.core)) if(nb.core == 1) parallel <- FALSE
  if(is.null(nb.core) & !parallel) nb.core <- 1
  
  if(method.abc != "rejection") {
    stop("only abc method rejection is implemented - ongoing work",
         call. = FALSE)
  }
  
  if(!is.null(pool)) {
    # Avoid factors in columns of pool
    sel <- which(lapply(pool, class) == "factor")
    for(i in sel) pool[,i] <- as.character(pool[,i])
  }
  
  if (ncol(pool) >= 3) {
    if (length(formals(f.sumstats))>1) {
      # In this case, traits needs to be defined
      # Compute mean or most frequent trait value per species
      sel <- which(lapply(pool, class) == "numeric") 
      traits <- data.frame(
        sapply(3:ncol(pool), function(x) 
          if(sel[x]) {
            tapply(data.frame(pool[,x]), pool[, 2], mean)
          } else tapply(data.frame(pool[,x]), pool[, 2], function(y) names(y[which.max(table(y))]))), 
        stringsAsFactors = F)
    }
  } else if (is.null(traits) & ncol(pool) < 3) 
    warning("Trait information is not provided", call. = FALSE)
  
  # Community size
  if (!(multi %in% c("single", "tab", "seqcom"))){
    stop("multi parameter must be either single, tab or seqcom.", call. = FALSE)
  }
  
  if (multi == "single") {
    J <- nrow(comm.obs)
    nb.com <- 1
  } else {
    if (multi == "tab") {
      # comm.obs is a species-by-site matrix/data.frame
      J <- apply(comm.obs, 1, function(x) sum(x, na.rm = TRUE))
      # if the dataset includes relative proportions, the columns must sum to 1
      if(prop & any(J != 1)) {
        stop("Relative species abundances must sum to 1", call. = FALSE)
      }
      nb.com <- nrow(comm.obs)
    } else if (multi == "seqcom") {
      # comm.obs is a list of communities with individuals on rows in each
      # community
      J <- lapply(comm.obs, nrow)
      nb.com <- length(comm.obs)
    }
  }
  
  if(multi == "seqcom"){
    if (length(formals(f.sumstats)) == 1) {
      stats.obs <- lapply(comm.obs, function(x) as.vector(f.sumstats(x)))
    } else {
      stats.obs <- lapply(comm.obs,
                          function(x) as.vector(f.sumstats(x,traits)))
    }
  } else {
    if (length(formals(f.sumstats)) == 1) {
      stats.obs <- as.vector(f.sumstats(comm.obs))
    } else {
      stats.obs <- as.vector(f.sumstats(comm.obs, traits))
    }
  }
  
  # Defining a list with the bounds of prior distributions
  prior=c();
  if(is.null(filt.abc)) {
    prior[[length(prior)+1]]  <- c("unif",0.1,1)
  }else{
    for(i in 1:nrow(par.filt)) prior[[i]] <- c("unif",par.filt[i,1],par.filt[i,2])
    prior[[length(prior)+1]] <- c("unif",0.1,1) # m 
    if(prop) 
    {
      prior[[length(prior)+1]] <- c("unif",100,1000)
      warning("The prior of community size is uniform between 100 and 1000",
              call. = FALSE)
    }
  }
  
  coalesc_model <- function(par, traits, prop, J, pool, filt.abc, f.sumstats, parallel) {
    stats.samp <- NA
    try({
      if(is.null(filt.abc)){
        if(!prop){
          comm.samp <- coalesc(J, m = par[length(par)], filt = NULL,
                               add = F,  var.add =NULL, pool=pool, 
                               traits = NULL,
                               checks = F)
          if (length(formals(f.sumstats))==1) {
            stats.samp <- as.vector(f.sumstats(comm.samp$com))
          } else {
            stats.samp <- as.vector(f.sumstats(comm.samp$com, traits))
          }
        }
      }else{
        if(!prop) {
          if(parallel){
            set.seed(par[1])
            comm.samp <- coalesc(J, m = par[length(par)], filt = function(x) filt.abc(x, par[2:(length(par)-1)]),
                                 add = F,  var.add = NULL, pool = pool, 
                                 traits = NULL, Jpool = 50 * J, verbose = FALSE, checks = F)
            if (length(formals(f.sumstats))==1) {
              stats.samp <- as.vector(f.sumstats(comm.samp$com))
            } else {
              stats.samp <- as.vector(f.sumstats(comm.samp$com, traits))
            }} else {
              comm.samp <- coalesc(J, m = par[length(par)], filt = function(x) filt.abc(x, par[1:(length(par)-1)]),
                                   add = F,  var.add = NULL, pool = pool, 
                                   traits = NULL, checks = F)
              if (length(formals(f.sumstats))==1) {
                stats.samp <- as.vector(f.sumstats(comm.samp$com))
              } else {
                stats.samp <- as.vector(f.sumstats(comm.samp$com, traits))
              }
            }
        } else {
          comm.samp <- coalesc(
            par[length(par)], m = par[length(par)-1], 
            filt = function(x) filt.abc(x, par[-((length(par)-1):length(par))]), 
            pool = pool, traits = traits, checks = F)
          comm.samp$com <- t(table(comm.samp$com[,2])/par[length(par)])
        }
        if (length(formals(f.sumstats)) == 1) {
          stats.samp <- as.vector(f.sumstats(comm.samp$com))
        } else {
          stats.samp <- as.vector(f.sumstats(comm.samp$com, traits))
        }
      }
      
    })
    return(stats.samp)
  }
  
  # For debug
  #test <- c()
  #for(i in 1:100) test[[i]] <- coalesc_model(c(runif(0,1),runif(0,1),runif(0,1)), traits, J, pool, filt.abc, f.sumstats)
  
  #if(multi!="single") stop("multi option is not implemented - ongoing work")
  
  set.seed(1)
  
  if(type == "seq")
  {
    if(method.seq == "Lenormand"){
      pacc = 0.05 # Can be set by user (to be included in input)
      res.abc <- EasyABC::ABC_sequential(
        method              = method.seq, 
        model               = function(x) {
          coalesc_model(x, traits,prop,J,pool,filt.abc,f.sumstats, parallel)
          },
        prior               = prior, 
        nb_simul            = nb.samp,
        summary_stat_target = stats.obs, 
        p_acc_min           = pacc, 
        use_seed            = FALSE, 
        n_cluster           = nb.core,
        alpha               = alpha)
    } else if(method.seq=="Beaumont") 
    { stop("Beaumont method not implemented - ongoing work", call. = FALSE)
      tol_tab <- c(tol,tol/2,tol/5)
      res.abc <- EasyABC::ABC_sequential(method=method.seq,
                                         model=function(par) coalesc_model(par, traits, prop, J, pool, filt.abc, f.sumstats, parallel),
                                         prior=prior,
                                         nb_simul=nb.samp, 
                                         summary_stat_target=stats.obs,
                                         tolerance_tab=tol_tab, 
                                         use_seed=F, 
                                         n_cluster=nb.core)
    } else if(method.seq=="Drovandi") 
    {stop("Drovandi method not implemented - ongoing work", call. = FALSE)
      res.abc <- EasyABC::ABC_sequential(method=method.seq, 
                                         model=function(par) coalesc_model(par, traits, prop, J, pool, filt.abc, f.sumstats, parallel),
                                         prior=prior, 
                                         nb_simul=nb.samp, summary_stat_target=stats.obs, first_tolerance_level_auto = T, use_seed=F, 
                                         n_cluster=nb.core)
    }
    
  } else if(type=="mcmc")
  { if(method.mcmc == "Marjoram_original"){
    res.abc <- EasyABC::ABC_mcmc(
      method=method.mcmc, 
      model=function(x) coalesc_model(x, traits, prop, J, pool,filt.abc, f.sumstats, parallel),
      prior=prior, 
      summary_stat_target=stats.obs,
      n_cluster=nb.core,
      n_rec = nb.samp)} else {
        stop("mcmc methods other than Marjoram_original are not implemented - ",
             "ongoing work", call. = FALSE)
      }
    
  } else if(type=="annealing")
  { stop("SABC is not implemented - ongoing work", call. = FALSE)
    #res <- EasyABC::SABC(r.model, r.prior, d.prior, n.sample, eps.init, iter.max,
    #     v=ifelse(method=="informative",0.4,1.2), beta=0.8,
    #     delta=0.1, resample=5*n.sample, verbose=n.sample,
    #     method="noninformative", adaptjump=TRUE,
    #     summarystats=FALSE, y=NULL, f.summarystats=NULL)
  }
  
  if(type != "standard"){
    comm.abc <- list(obs = stats.obs,
                     abc = res.abc,
                     call = match.call()) }
  return(comm.abc)
}

coalesc_abc_std <- function(comm.obs, pool = NULL, multi = "single", prop = F, traits = NULL,
                            f.sumstats, filt.abc = NULL, migr.abc = NULL, size.abc = NULL, add = F,
                            var.add = NULL, params = NULL, par.filt = NULL, par.migr = NULL, 
                            par.size = NULL, constr = NULL, scale = F, dim.pca = NULL, svd = F, theta.max = NULL, 
                            nb.samp = 10^6, parallel = TRUE, nb.core = NULL, tol = NULL, pkg = NULL, 
                            method.abc = "rejection")
{
  initial_checks(comm.obs = comm.obs, pool = pool, multi = multi, prop = prop,
                 traits = traits, f.sumstats = f.sumstats, filt.abc = filt.abc,
                 migr.abc = migr.abc, size.abc = size.abc,  add = add,
                 var.add = var.add, params = params, par.filt = par.filt,
                 par.migr = par.migr, par.size = par.size, constr = constr,
                 scale = scale, dim.pca = dim.pca, svd = svd,
                 theta.max = theta.max, nb.samp = nb.samp, parallel = parallel,
                 nb.core = nb.core, tol = tol, pkg = pkg,
                 method.abc = method.abc)
  
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
  sim <- do.simul.coalesc(J, pool, multi, prop, nb.com, traits, f.sumstats, filt.abc, migr.abc,
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
                             traits = NULL, f.sumstats = NULL, filt.abc = NULL, migr.abc = NULL, 
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
                          theta = params.samp[length(params.samp)], checks = F)$com
        } else {
          pool <- coalesc(mean(J)*100,
                          theta = params.samp[length(params.samp)], checks = F)$com
        }
        params.samp <- params.samp[-length(params.samp)]
      }
      
      ifelse(!is.null(filt.abc),
             filt <- ifelse(!add, 
                            function(x) filt.abc(x, params.samp[1:nrow(par.filt)]),
                            function(x, var.add) filt.abc(x, params.samp[1:nrow(par.filt)], var.add)),
             filt <- NULL)
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
                                                    pool = pool.loc, traits = traits, checks = F)$com
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
                                         pool = pool, traits = traits,
                                         checks = F)
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
        stop("Package ecolottery is not available", call. = FALSE)
      }
      # Other required packages
      if (!is.null(pkg)) for (i in 1:length(pkg)) {
        if (!requireNamespace(pkg[i], quietly = TRUE)) {
          stop(paste("Package ", pkg[i], " is not available", sep = ""),
               call. = FALSE)
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

generate_prior <- function(pool = NULL, prop = F, constr = NULL, 
                           params = NULL, par.filt = NULL, par.migr = NULL,
                           par.size = NULL, theta.max = NULL, nb.samp = 10^6)
{
  # Function defining constraints on parameter values
  constr.test <- function(par.names, par.val) 
  {
    res <- T
    if(!is.null(constr)) {
      for(j in 1:length(par.names))
        assign(par.names[j], par.val[j], envir = environment())
      for(j in 1:length(constr)) 
        res <- res & lazyeval::lazy_eval(constr[j], environment())
    }
    return(res)
  }
  
  if (!is.null(constr) & !requireNamespace("lazyeval", quietly = TRUE)) {
    stop("coalesc_abc requires package lazy_eval to be installed",
         call. = FALSE)
  }
  
  # Uniform prior distributions of parameters
  prior <- list()
  dim.prior <- ifelse(!is.null(nrow(par.filt)), nrow(par.filt), 0) + 
    max(ifelse(!is.null(par.migr), nrow(par.migr), 0), 1) + 
    prop*max(ifelse(!is.null(par.size), nrow(par.size), 0), 1) + 
    as.numeric(is.null(pool))
  length(prior) <- dim.prior
  stop <- 0
  while(length(prior[[1]]) < nb.samp & stop < 100) {
    samp <- nb.samp - length(prior[[1]])
    i <- 0
    if(!is.null(par.filt)) {
      for (i in 1:nrow(par.filt)) {
        prior[[i]] <- c(prior[[i]], runif(samp, min = par.filt[i, 1],
                                          max = par.filt[i, 2]))
      }
    }
    if(!is.null(par.migr)) {
      for (j in (i+1):(i+nrow(par.migr))) {
        prior[[j]] <- c(prior[[j]], runif(samp, min = par.migr[j-i, 1], 
                                          max = par.migr[j-i, 2]))
        i <- j
      }
    } else {
      prior[[i+1]] <- c(prior[[i+1]], runif(samp, min = 0, max = 1))
      i <- i + 1
      j <- i
    }
    if(!is.null(par.migr)) {
      names(prior) <- c(rownames(par.filt), rownames(par.migr))
    } else {
      names(prior) <- c(rownames(par.filt), "m")
    }
    
    # TODO - Defining a lower bound at 0 for m and theta can entail issues when
    # species richness is 1 in simulated community; we should allow the user to
    # define the prior for m and theta in the future
    
    if(prop) 
    {
      if(!is.null(par.size))
      {
        for (i in (j+1):(j + nrow(par.size))) {
          prior[[i]] <- c(prior[[i]], runif(samp, min = par.size[i-j, 1],
                                            max = par.size[i-j, 2]))
        }
        names(prior)[(j+1):(j+nrow(par.size))] <- rownames(par.size)
      } else {
        prior[[j+1]] <- c(prior[[j+1]], runif(samp, min = 100, max = 1000))
        i <- j+1
        warning("The prior of community size is uniform between 100 and 1000",
                call. = FALSE)
        names(prior)[i] <- "J"
      }
    }
    
    if (is.null(pool)) {
      prior[[i + 1]] <- c(prior[[i + 1]], runif(samp, min = 0, max = theta.max))
      names(prior)[i+1] <- "theta"
    }
    
    constr.sel <- sapply(1:nb.samp, function(x) {
      constr.test(names(prior), unlist(lapply(prior, function(y) y[x])))
    })
    prior <- lapply(prior, function(x) x[constr.sel])
    
    stop <- stop + 1
  }
  
  if(stop==100) warning("Fail to meet constraints for nb.samp sets of ",
                        "parameter values", call. = FALSE)
  
  return(prior)
}

initial_checks <- function(comm.obs = NULL, pool = NULL, multi = "single", prop = F, traits = NULL, 
                           f.sumstats = NULL, filt.abc = NULL,  migr.abc = NULL, size.abc = NULL, 
                           add = NULL, var.add = NULL, params = NULL, par.filt = NULL, par.migr = NULL, 
                           par.size = NULL, constr = NULL, scale = F, dim.pca = NULL, svd = NULL,
                           theta.max = NULL, nb.samp = 10^6, parallel = F, nb.core = NULL, tol = NULL, type = "standard",
                           method.seq = "Lenormand", method.mcmc = "Marjoram_original",
                           method.abc = "rejection", alpha = 0.5, pkg = NULL)
{
  # Checks on f.sumstats
  if(!is.function(f.sumstats)) {
    stop("You must provide a function to calculate summary statistics",
         "(f.sumstats)", call. = FALSE)
  }
  if (length(formals(f.sumstats)) > 3) {
    stop("f.sumstats must be a function of up to three arguments", call. = FALSE)
  }
  
  if(is.null(filt.abc)){
    warning("No habitat filtering function provided. Neutral communities ",
            "will be simulated and only m will be estimated", call. = FALSE)
  }
  
  # Should ABC analysis be performed
  if (is.null(tol)){
    warning("You must provide a tolerance value for ABC computation.\n",
            "The function will only provide simulations and will not perform ",
            "ABC analysis.", call. = FALSE)
  } 
  else {
    if (!requireNamespace("abc", quietly = TRUE)) {
      stop("coalesc_abc requires package abc to be installed", call. = FALSE)
    }
    if(!is.null(method.abc)) if(!method.abc%in%c("rejection", "loclinear", "neuralnet", "ridge"))
      stop("method.abc should be either rejection, loclinear, neuralnet ",
           "or ridge", call. = FALSE)
  }
  
  # Checks on multi
  if (!(multi %in% c("single", "tab", "seqcom"))){
    stop("multi parameter must be either single, tab or seqcom.", call. = FALSE)
  }
  if (!multi == "single" & class(pool) == "list") {
    if((multi == "tab" & length(pool) != nrow(comm.obs)) | (multi == "seqcom" & length(pool) != length(comm.obs)))
      stop("When several pools are given (list), there must be one for each ",
           "community", call. = FALSE)
  }
  if (multi=="single" & class(pool)=="list" & length(pool)>1)
    stop("If multi=single, pool should not be a list of more than one element",
         call. = FALSE)
  
  if (multi=="tab") if (any(rowSums(comm.obs)==0)) {
    stop("There should not be communities with 0 individuals", call. = FALSE)
  }
  
  if(prop & multi!="tab") {
    stop("prop data can only be handled in tab format", call. = FALSE)
  }
  
  if (is.null(traits)){
    warning("Trait information is not provided", call. = FALSE)
  }
  
  # Other required packages
  for (i in pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package ", pkg, " is not available", sep = ""), call. = FALSE)
    }   
  }
}

