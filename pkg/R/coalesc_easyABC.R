coalesc_easyABC <- function(comm.obs, pool = NULL, multi = "single", prop = F, traits = NULL,
                            f.sumstats, filt.abc = NULL, filt.vect = F, migr.abc = NULL, size.abc = NULL, 
                            add = F, var.add = NULL, params = NULL, par.filt = NULL, par.migr = NULL, 
                            par.size = NULL, constr = NULL, scale = F, dim.pca = NULL, svd = F, theta.max = NULL, 
                            nb.samp = 10^6, parallel = TRUE, nb.core = NULL, tol = NULL, pkg = NULL, 
                            type = NULL, method.seq = "Lenormand", method.mcmc = "Marjoram_original", method.abc = "rejection",
                            alpha = 0.5) 
{
  if(!is.null(migr.abc) | !is.null(size.abc)) {
    migr.abc <- NULL
    size.abc <- NULL
    warning("migr.abc, size.abc not available with other sampling approach ",
            "than standard", call. = FALSE)
  }
  if(!is.null(constr)) {
    warning("constr not available with other sampling approach ",
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
  }
  
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
  nb.sumstats <- length(stats.obs)
  
  # Defining a list with the bounds of prior distributions
  # By default, a uniform prior is defined for migration rate, between 0 and 1
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
  
  coalesc_model <- function(par, traits, prop, J, pool, filt.abc, filt.vect, f.sumstats, 
                            parallel) {
    stats.samp <- NA
    try({
      if(is.null(filt.abc)){
        if(!prop){
          comm.samp <- coalesc(J, m = par[length(par)], filt = NULL,
                               add = F,  var.add = NULL, pool = pool, 
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
                                 filt.vect = filt.vect,
                                 add = F,  var.add = NULL, pool = pool, 
                                 traits = NULL, Jpool = 50 * J, verbose = FALSE, checks = F)
            if (length(formals(f.sumstats))==1) {
              stats.samp <- as.vector(f.sumstats(comm.samp$com))
            } else {
              stats.samp <- as.vector(f.sumstats(comm.samp$com, traits))
            }} else {
              comm.samp <- coalesc(J, m = par[length(par)], filt = function(x) filt.abc(x, par[1:(length(par)-1)]),
                                   filt.vect = filt.vect,
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
            filt.vect = filt.vect,
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
  #for(i in 1:100) test[[i]] <- coalesc_model(c(runif(0,1),runif(0,1),runif(0,1)), 
  #traits, J, pool, filt.abc, filt.vect, f.sumstats)
  
  #if(multi!="single") stop("multi option is not implemented - ongoing work")
  
  set.seed(1)
  
  if(type == "seq")
  {
    if(method.seq == "Lenormand"){
      pacc = 0.05 # Can be set by user (to be included in input)
      res.abc <- EasyABC::ABC_sequential(
        method              = method.seq, 
        model               = function(x) {
          coalesc_model(x, traits, prop, J, pool, filt.abc, 
                        filt.vect, f.sumstats, parallel)
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
                                         model=function(par) coalesc_model(par, traits, prop, 
                                                                           J, pool, filt.abc, 
                                                                           filt.vect, f.sumstats, parallel),
                                         prior=prior,
                                         nb_simul=nb.samp, 
                                         summary_stat_target=stats.obs,
                                         tolerance_tab=tol_tab, 
                                         use_seed=F, 
                                         n_cluster=nb.core)
    } else if(method.seq=="Drovandi") 
    {stop("Drovandi method not implemented - ongoing work", call. = FALSE)
      res.abc <- EasyABC::ABC_sequential(method=method.seq, 
                                         model=function(par) coalesc_model(par, traits, prop, J, pool, filt.abc, filt.vect, f.sumstats, parallel),
                                         prior=prior, 
                                         nb_simul=nb.samp, summary_stat_target=stats.obs, first_tolerance_level_auto = T, use_seed=F, 
                                         n_cluster=nb.core)
    }
    
  } else if(type=="mcmc")
  { if(method.mcmc == "Marjoram_original"){
    res.abc <- EasyABC::ABC_mcmc(
      method=method.mcmc, 
      model=function(x) coalesc_model(x, traits, prop, J, pool,filt.abc, filt.vect, f.sumstats, parallel),
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
