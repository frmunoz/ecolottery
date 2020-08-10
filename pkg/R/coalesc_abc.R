coalesc_abc <- function(comm.obs, pool = NULL, multi = "single", prop = FALSE,
                        traits = NULL, f.sumstats, filt.abc = NULL, filt.vect = F,
                        migr.abc = NULL, size.abc = NULL, add = FALSE,
                        var.add = NULL, params  = NULL, par.filt = NULL,
                        par.migr = NULL, par.size = NULL, constr = NULL,
                        scale = FALSE, dim.pca = NULL, svd = FALSE, theta.max = NULL,
                        nb.samp = 10^6, parallel = FALSE, nb.core = NULL,
                        tol = NULL, type = "standard", method.seq = "Lenormand",
                        method.mcmc = "Marjoram_original",
                        method.abc = "rejection", alpha = 0.5, pkg = NULL) 
{
  # Work in global environment
  globalenv()
  
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
  
  initial_checks(comm.obs = comm.obs, pool = pool, multi = multi, prop = prop,
                 traits = traits, f.sumstats = f.sumstats,
                 filt.abc = filt.abc, filt.vect = filt.vect, migr.abc = migr.abc,
                 size.abc = size.abc, params = params, par.filt = par.filt,
                 par.migr = par.migr,  par.size = par.size, scale = scale,
                 theta.max = theta.max, nb.samp = nb.samp,
                 parallel = parallel, nb.core = nb.core, tol = tol,
                 type = type, method.seq = method.seq,
                 method.mcmc = method.mcmc, method.abc = method.abc,
                 alpha = alpha, pkg = pkg)
  
  # The options migr.abc = NULL, size.abc = NULL, add = F, 
  # var.add = NULL, par.filt = NULL, par.migr = NULL, 
  # par.size = NULL, constr = NULL, scale = T, dim.pca = NULL, svd = F, pkg = NULL
  # need to be included
  if (type == "standard")
    return(coalesc_abc_std(comm.obs, pool, multi, prop, traits, f.sumstats,
                           filt.abc, filt.vect, migr.abc, size.abc, add, var.add, params,
                           par.filt, par.migr, par.size, constr,
                           scale, dim.pca, svd, theta.max, nb.samp,
                           parallel, nb.core, tol, pkg, method.abc)) 
  else return(coalesc_easyABC(comm.obs, pool, multi, prop, traits, f.sumstats,
                              filt.abc, filt.vect, migr.abc, size.abc, add, var.add, params,
                              par.filt, par.migr, par.size, constr,
                              scale, dim.pca, svd, theta.max, nb.samp,
                              parallel, nb.core, tol, pkg, type, method.seq, 
                              method.mcmc, method.abc, alpha)) 
}

initial_checks <- function(comm.obs = NULL, pool = NULL, multi = "single", prop = F, traits = NULL, 
                           f.sumstats = NULL, filt.abc = NULL, filt.vect = F, migr.abc = NULL, size.abc = NULL, 
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

  # Other required packages
  for (i in pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package ", pkg, " is not available", sep = ""), call. = FALSE)
    }   
  }
}
