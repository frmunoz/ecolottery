coalesc_abc2 <- function (comm.obs, pool, multi = "single", prop = F, traits = NULL, f.sumstats, filt.abc, params,
                          theta.max = NULL, nb.samp = 10^6, parallel = T, nb.core = NULL, tol = NULL, type = "standard", 
                          method.seq = "Lenormand", method.mcmc = "Marjoram_original", method.abc = NULL, scale = F) 
{
  # This alternative function uses the sequential algorithms provided in EasyABC
  if(is.null(type)) error("Standard ABC analysis with coalesc_abc")
  
  if (is.null(pool)) {
    warning(paste0("No species pool provided: pool will be simulated ",
                   "and logseries theta estimated"))
    if (is.null(theta.max))
      theta.max <- 500
    stop("Estimation of theta not implemented - Ongoing work")
  }
  
  if(prop & multi!="tab")
    stop("prop data can only be handled in tab format")
  
  if (!requireNamespace("abc", quietly = TRUE)) 
    stop("coalesc_abc requires package abc to be installed")
  
  if (!requireNamespace("EasyABC", quietly = TRUE)) 
    stop("coalesc_abc requires package EasyABC to be installed")
  #for (i in 1:length(pkg)) if (!requireNamespace(pkg[i], quietly = TRUE)) 
  #  stop(paste("Package ", pkg[i], " is not available", sep = ""))

  if(!type%in%c("seq","mcmc","annealing","standard"))
    stop("type should be either standard, seq, mcmc or annealing")
  
  if(type=="seq") if(!method.seq%in%c("Beaumont", "Drovandi", "Delmoral", "Lenormand", "Emulation"))
    stop("method.seq should be either Beaumont, Drovandi, Delmoral, Lenormand or Emulation")
  
  if(type=="mcmc") if(!method.mcmc%in%c("Marjoram_original", "Marjoram", "Wegmann"))
    stop("method.mcmc should be either Marjoram_original, Marjoram or Wegmann")
  
  if(!is.null(nb.core) & nb.core == 1) parallel <- F
  
  if(!is.null(method.abc)) if(!method.abc%in%c("rejection", "loclinear", "neuralnet", "ridge"))
    stop("method.abc should be either rejection, loclinear, neuralnet or ridge")
  
  if (ncol(pool) >= 3) 
    traits <- data.frame(apply(data.frame(pool[, -(1:2)]), 2, function(x) tapply(x, pool[, 2], mean))) else if (is.null(traits) & ncol(pool) < 3) 
    warning("Trait information is not provided")
  
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
  
  if(multi == "seqcom"){
    if (length(formals(f.sumstats))==1) {
      stats.obs <- lapply(comm.obs, function(x) as.vector(f.sumstats(x)))
    } else {
      stats.obs <- lapply(comm.obs, function(x) as.vector(f.sumstats(x, traits)))
    }
  } else {
    if (length(formals(f.sumstats))==1) {
      stats.obs <- as.vector(f.sumstats(comm.obs))
    } else {
      stats.obs <- as.vector(f.sumstats(comm.obs, traits))
    }
  }
  
  # Definition of prior distribution to be improved
  prior=c();
  for(i in 1:nrow(params)) prior[[i]] <- c("unif",params[i,1],params[i,2])
  prior[[length(prior)+1]] <- c("unif",0.1,1)
  if(prop) 
  {
    prior[[length(prior)+1]] <- c("unif",100,1000)
    warning("The prior of community size is uniform between 100 and 1000")
  }

  coalesc_model <- function(par, traits, prop, J, pool, filt.abc, f.sumstats) {
      stats.samp <- NA
      try({
            if(!prop) {
              comm.samp <- coalesc(J, m = par[length(par)], 
                                 filt = function(x) filt.abc(x, par[-length(par)]), 
                                 pool = pool, traits = traits)
            } else {
              comm.samp <- coalesc(par[length(par)], m = par[length(par)-1], 
                                      filt = function(x) filt.abc(x, par[-((length(par)-1):length(par))]), 
                                      pool = pool, traits = traits)
              comm.samp$com <- t(table(comm.samp$com[,2])/par[length(par)])
            }
            if (length(formals(f.sumstats))==1) {
              stats.samp <- as.vector(f.sumstats(comm.samp$com))
            } else {
              stats.samp <- as.vector(f.sumstats(comm.samp$com, traits))
            }
          })
      return(stats.samp)
  }
  
  # For debug
  #test <- c()
  #for(i in 1:100) test[[i]] <- coalesc_model(c(runif(0,1),runif(0,1),runif(0,1)), traits, J, pool, filt.abc, f.sumstats)
  
  #if(multi!="single") stop("multi option is not implemented - ongoing work")
  
  set.seed(1)
  if(type == "standard"){ # performing standard abc with chosen abc.method 
    comm.sim <- ABC_rejection(model = function(x) coalesc_model(x, traits, prop,J,pool,filt.abc,f.sumstats),
                              prior = prior, 
                              nb_simul = nb.samp, 
                              n_cluster = nb.core)
    if(scale){
      
      stats.mean <- apply(comm.sim$stats,2, function(x) mean(x, na.rm =T))
      stats.sd <- apply(comm.sim$stats,2, function(x) sd(x, na.rm =T))
      comm.sim$stats.scaled <- t(apply(comm.sim$stats, 1, function(x) (x - stats.mean/stats.sd)))
      stats.obs.scaled <- (stats.obs - stats.mean) / stats.sd
      
      res.abc <- abc::abc(stats.obs.scaled, comm.sim$param, comm.sim$stats.scaled, tol= tol, method = method.abc)
      
      comm.abc <- list(par = comm.sim$param,
                       obs = stats.obs,
                       obs.scaled = stats.obs.scaled,
                       ss = comm.sim$stats,
                       ss.scaled = comm.sim$stats.scaled,
                       scale = t(data.frame(mean = stats.mean, sd = stats.sd)),
                       abc = res.abc)
      
    } else {res.abc <- abc::abc(stats.obs, comm.sim$param, comm.sim$stats, tol= tol, method = method.abc) 
    
    comm.abc <- list(par = comm.sim$param,
                     obs = stats.obs,
                     ss = comm.sim$stats,
                     abc = res.abc) }
  }
  
  if(type=="seq")
  {
    if(method.seq=="Lenormand"){
      pacc=0.05 # Can be set by user (to be included in input)
      res.abc <- EasyABC::ABC_sequential(method=method.seq, 
                                     model=function(par) coalesc_model(par, traits, prop, J, pool, filt.abc, f.sumstats),
                                     prior=prior, 
                                     nb_simul=nb.samp,
                                     summary_stat_target=stats.obs, 
                                     p_acc_min=pacc, 
                                     use_seed=F, 
                                     n_cluster=nb.core)
    } else if(method.seq=="Beaumont") 
      { stop("Beaumont method not implemented - ongoing work")
      tol_tab <- c(tol,tol/2,tol/5)
      res.abc <- EasyABC::ABC_sequential(method=method.seq,
                                         model=function(par) coalesc_model(par, traits, prop, J, pool, filt.abc, f.sumstats),
                                         prior=prior,
                                         nb_simul=nb.samp, 
                                         summary_stat_target=stats.obs,
                                         tolerance_tab=tol_tab, 
                                         use_seed=F, 
                                         n_cluster=nb.core)
    } else if(method.seq=="Drovandi") 
    {stop("Drovandi method not implemented - ongoing work")
     res.abc <- EasyABC::ABC_sequential(method=method.seq, 
                                    model=function(par) coalesc_model(par, traits, prop, J, pool, filt.abc, f.sumstats),
                                    prior=prior, 
                                    nb_simul=nb.samp, summary_stat_target=stats.obs, first_tolerance_level_auto = T, use_seed=F, 
                                    n_cluster=nb.core)
    }
    
  } else if(type=="mcmc")
  { if(method.mcmc == "Marjoram_original"){
    res.abc <- EasyABC::ABC_mcmc(method=method.mcmc, 
                                 model=function(x) coalesc_model(x, traits, prop, J, pool,filt.abc, f.sumstats),
                                 prior=prior, 
                                 summary_stat_target=stats.obs,
                                 n_cluster=nb.core,
                                 n_rec = nb.samp) } else stop("mcmc methods other than Marjoram_original are not implemented - ongoing work")
    
  } else if(type=="annealing")
  { stop("SABC is not implemented - ongoing work")
    #res <- EasyABC::SABC(r.model, r.prior, d.prior, n.sample, eps.init, iter.max,
    #     v=ifelse(method=="informative",0.4,1.2), beta=0.8,
    #     delta=0.1, resample=5*n.sample, verbose=n.sample,
    #     method="noninformative", adaptjump=TRUE,
    #     summarystats=FALSE, y=NULL, f.summarystats=NULL)
  }
  
  if(type != "standard"){
    comm.abc <- list(obs = stats.obs,
                   abc = res.abc) }
  return(comm.abc)
}
