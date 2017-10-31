coalesc_abc2 <- function (comm.obs, pool, multi = "single", traits = NULL, f.sumstats, filt.abc, params,
                          theta.max = NULL, nb.samp = 10^6, parallel = T, tol = NULL, type = "seq", 
                          method.seq = "Beaumont", method.mcmc = "Marjoram", method.abc = NULL) 
{
  # This alternative function uses the sequential algorithms provided in EasyABC
  
   if (is.null(pool)) {
    warning(paste0("No species pool provided: pool will be simulated ",
                   "and logseries theta estimated"))
    if (is.null(theta.max))
      theta.max <- 500
    stop("Estimation of theta not implemented - Ongoing work")
  }
  
  if (!requireNamespace("abc", quietly = TRUE)) 
    stop("coalesc_abc requires package abc to be installed")
  
  if (!requireNamespace("EasyABC", quietly = TRUE)) 
    stop("coalesc_abc requires package EasyABC to be installed")
  #for (i in 1:length(pkg)) if (!requireNamespace(pkg[i], quietly = TRUE)) 
  #  stop(paste("Package ", pkg[i], " is not available", sep = ""))

  if(!type%in%c("seq","mcmc","annealing"))
    stop("type should be either seq, mcmc or annealing")
  
  if(type=="seq") if(!method.seq%in%c("Beaumont", "Drovandi", "Delmoral", "Lenormand", "Emulation"))
    stop("method.seq should be either Beaumont, Drovandi, Delmoral, Lenormand or Emulation")
  
  if(type=="mcmc") if(!method.mcmc%in%c("Marjoram_original", "Marjoram", "Wegmann"))
    stop("method.mcmc should be either Marjoram_original, Marjoram or Wegmann")
  
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
  
  coalesc_model <- function(par, traits, J, pool, filt.abc, f.sumstats) {
      stats.samp <- NA
      try({
            comm.samp <- coalesc(J, m = par[length(par)], 
                                 filt = function(x) filt.abc(x, par[-length(par)]), 
                                 pool = pool, traits = traits)
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
  
  if(multi!="single") stop("multi option is not implemented - ongoing work")
  
  set.seed(1)

  if(type=="seq")
  {
    if(method.seq=="Lenormand") 
    {
      pacc=0.05 # Can be set by user (to be included in input)
      res <- EasyABC::ABC_sequential(method=method.seq, model=function(par) coalesc_model(par, traits, J, pool, 
          filt.abc, f.sumstats), prior=prior, nb_simul=nb.samp, summary_stat_target=stats.obs, 
          p_acc_min=pacc, use_seed=F)
    } else if(method.seq=="Beaumont") 
      {
      tol_tab <- c(tol,tol/2,tol/5)
      res <- EasyABC::ABC_sequential(method=method.seq, 
          model=function(par) coalesc_model(par, traits, J, pool, filt.abc, f.sumstats), prior=prior, 
          nb_simul=nb.samp, summary_stat_target=stats.obs, tolerance_tab=tol_tab, use_seed=F)
    } else if(method.seq=="Drovandi") 
    {
     res <- EasyABC::ABC_sequential(method=method.seq, 
                                     model=function(par) coalesc_model(par, traits, J, pool, filt.abc, f.sumstats), prior=prior, 
                                     nb_simul=nb.samp, summary_stat_target=stats.obs, first_tolerance_level_auto = T, use_seed=F)
    }
    
  } else if(type=="mcmc")
  {
    res <- EasyABC::ABC_mcmc(method=method.mcmc, model=function(par) coalesc_model(par, traits, J, pool, 
          filt.abc, f.sumstats), prior=prior, summary_stat_target=stats.obs)
  } else if(type=="annealing")
  {
    stop("SABC is not implemented - ongoing work")
    #res <- EasyABC::SABC(r.model, r.prior, d.prior, n.sample, eps.init, iter.max,
    #     v=ifelse(method=="informative",0.4,1.2), beta=0.8,
    #     delta=0.1, resample=5*n.sample, verbose=n.sample,
    #     method="noninformative", adaptjump=TRUE,
    #     summarystats=FALSE, y=NULL, f.summarystats=NULL)
  }
  
  if(!is.null(method.abc))
  {
    res.abc <- abc(stats.obs, res$param, res$stats, tol=tol, method=method.abc)
    return(list(res=res, stats.obs=stats.obs, abc=res.abc))
  } else return(list(res=res, stats.obs=stats.obs))
}
