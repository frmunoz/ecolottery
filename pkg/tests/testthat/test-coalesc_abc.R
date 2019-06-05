context("Test coalesc_abc() structure")

test_that("coalesc_abc() works and returns desired result structure", {
  set.seed(1)
  
  # 1/ Analysis of community assembly in a single community
  
  # Trait-dependent filtering function
  filt_gaussian <- function(t, par) exp(-(t - par[1])^2/(2*par[2]^2))
  
  # Definition of parameters of this function and their range
  par.range <- data.frame(rbind(c(0, 1), c(0.05, 1)))
  row.names(par.range) <- c("topt", "sigmaopt")
  
  # Basic summary statistics describing local community composition
  f.sumstats <- function(com) array(dimnames=list(c("cwm", "cwv", "cws",
                                                    "cwk", "S", "Es")),
                                    c(mean(com[,3]), var(com[,3]), 
                                      e1071::skewness(com[,3]),  
                                      e1071::kurtosis(com[,3]),
                                      vegan::specnumber(table(com[,2])),
                                      vegan::diversity(table(com[,2]))))
  
  # An observed community is simulated with known parameter values
  comm <- coalesc(J = 400, m = 0.5, theta = 50,
                  filt = function(x) filt_gaussian(x, c(0.2, 0.1)))
 
  # ABC estimation of the parameters of filtering and migration based on observed 
  # community composition 
  # Number of values to sample in prior distributions
  nb.samp <- 3
  sink(tempfile())
  res.single <- coalesc_abc(comm$com, comm$pool, f.sumstats = f.sumstats,
                           filt.abc = filt_gaussian, par.filt = par.range, 
                           nb.samp = nb.samp, parallel = TRUE, tol = 0.5,
                           pkg = NULL, #c("e1071","vegan"),
                           method.abc = "neuralnet")
  sink()
  expect_is(res.single, "list")
  expect_named(res.single, c("par", "obs", "obs.scaled", "ss", "ss.scale", "abc", "call"))
  expect_is(res.single$par, "matrix")
  expect_is(res.single$obs, "array")
  expect_is(res.single$obs.scaled, "array")
  expect_is(res.single$ss, "matrix")
  expect_is(res.single$ss.scale, "data.frame")
  expect_is(res.single$abc, "abc")
  
  # Alternative option using a MCMC sampling approach
  sink(tempfile())
  res.single2 <- coalesc_abc(comm$com, comm$pool, f.sumstats=f.sumstats, filt.abc=filt_gaussian, 
                             par.filt=par.range, nb.samp=nb.samp, tol=0.1, type = "mcmc")
  sink()
  expect_is(res.single2, "list")
  expect_named(res.single2, c("obs", "abc", "call"))
  expect_is(res.single2$obs, "numeric")
  expect_is(res.single2$abc, "list")
  
  # 2/ Analysis of community assembly in multiple communities with a common pool of 
  # immigrants
  
  # When the input is a site-species matrix, use argument multi="tab"
  
  # 2.1/ Environment-dependent environmental filtering
  
  # The variation of optimal trait values can depend on environmental variation across communities
  filt_gaussian_env <- function(t, par, env) exp(-(t-par[1]*env+par[2])^2/(2*par[3]^2))
  # Vector of environmental conditions across communities
  env1 <- matrix(runif(20))
  
  # Simulate a set of 20 communities with varying environmental filtering in
  # different environmental conditions
  # A common source pool of immigrants
  pool <- coalesc(J = 10000, theta = 50)$com
  meta1 <- c()
  for(i in 1:20)
    meta1[[i]] <- coalesc(J = 400, pool = pool, m = 0.5, 
                          filt = function(x) filt_gaussian_env(x, c(0.5, 0.2, 0.1), env1[i,]))
  # Build a species-by-site table
  tab1 <- array(0, c(20, max(pool$sp)))
  for(i in 1:20) {
    compo <- table(meta1[[i]]$com$sp)
    tab1[i, as.numeric(names(compo))] <- compo
  }
  colnames(tab1) <- as.character(1:ncol(tab1))
  tab1 <- tab1[, colSums(tab1)!=0]
  
  # Definition of parameters and their range
  # We will estimate the slope a and intercept b of the linear relationship between 
  # optimal trait values and environmental conditions
  par.range <- data.frame(rbind(c(-1, 1), c(0, 1), c(0.05, 1)))
  row.names(par.range) <- c("a", "b", "sigmaopt")
  
  # Basic summary statistics
  # The function will provide trait-based statistics and taxonomic diversity in
  # each community
  f.sumstats <- function(tab, traits) array(dimnames=
                                              list(c(paste(rep("cwm",nrow(tab)),1:nrow(tab)), paste(rep("cwv",nrow(tab)),1:nrow(tab)), 
                                                     paste(rep("cws",nrow(tab)),1:nrow(tab)), paste(rep("cwk",nrow(tab)),1:nrow(tab)), 
                                                     paste(rep("S",nrow(tab)),1:nrow(tab)), paste(rep("Es",nrow(tab)),1:nrow(tab)), 
                                                     paste(rep("beta",nrow(tab)),1:nrow(tab)))),
                                            c(apply(tab, 1, function(x) Weighted.Desc.Stat::w.mean(traits[colnames(tab),], x)), 
                                              apply(tab, 1, function(x) Weighted.Desc.Stat::w.var(traits[colnames(tab),], x)), 
                                              apply(tab, 1, function(x) Weighted.Desc.Stat::w.skewness(traits[colnames(tab),], x)),  
                                              apply(tab, 1, function(x) Weighted.Desc.Stat::w.kurtosis(traits[colnames(tab),], x)),
                                              apply(tab, 1, function(x) sum(x!=0)),
                                              apply(tab, 1, vegan::diversity),
                                              colMeans(as.matrix(vegan::betadiver(tab,"w")))))
  
  # ABC estimation of the parameters of trait-environment relationship and of migration,
  # based on observed community composition 
  sink(tempfile())
  res.tab1 <- coalesc_abc(tab1, pool, multi = "tab", f.sumstats = f.sumstats,
                          filt.abc = filt_gaussian_env, par.filt = par.range, 
                          add = T, var.add = env1,
                          nb.samp = nb.samp, parallel = TRUE, tol = 0.5,
                          pkg = NULL,# c("e1071","vegan"),
                          method.abc = "neuralnet")
  sink()
  expect_is(res.tab1, "list")
  expect_named(res.tab1, c("par", "obs", "obs.scaled", "ss", "ss.scale", "abc", "call"))
  expect_is(res.tab1$par, "matrix")
  expect_is(res.tab1$obs, "array")
  expect_is(res.tab1$obs.scaled, "array")
  expect_is(res.tab1$ss, "matrix")
  expect_is(res.tab1$ss.scale, "data.frame")
  expect_is(res.tab1$abc, "abc")
  
  # 2.2/ Environment-dependent environmental filtering and immigration
  
  # In this case, the migration rate depends of another environmental variable
  # representing, e.g., community isolation
  # There will be two environmental variables used for parameter inference, one conditioning
  # environmental filtering, and the other conditioning migration
  env2 <-  matrix(runif(20))
  env <- cbind(env1, env2)
  colnames(env) <- c("env1", "env2")
  # Migration depends on environmental conditions following some linear relationship
  migr_env <- function(par, env) (par[2]-par[1])*env+par[1]
  
  # Simulate a set of 20 communities with environment-dependent environmental filtering 
  # and immigration
  meta2 <- c()
  for(i in 1:20)
    meta2[[i]] <- coalesc(J = 400, pool = pool, m = migr_env(c(0.25,0.75), env[i,2]), 
                          filt = function(x) filt_gaussian_env(x, c(0.5, 0.2, 0.1), env[i,1]))
  # Build a species-by-site table
  tab2 <- array(0, c(20, max(pool$sp)))
  for(i in 1:20) {
    compo <- table(meta2[[i]]$com$sp)
    tab2[i, as.numeric(names(compo))] <- compo
  }
  colnames(tab2) <- as.character(1:ncol(tab2))
  tab2 <- tab2[, colSums(tab2)!=0]
  
  # Definition of parameters and their range
  # We will estimate the slope a and intercept b of the linear relationship between 
  # optimal trait values and environmental variable 1, and slope c and intercept d in the
  # logistic function determining the variation of migration rate with environmental variable 2
  par.filt.range <- data.frame(rbind(c(-1, 1), c(0, 1), c(0.05, 1)))
  row.names(par.filt.range) <- c("a", "b", "sigmaopt")
  par.migr.range <- data.frame(rbind(c(0, 1), c(0, 1)))
  row.names(par.migr.range) <- c("c", "d")
  
  # ABC estimation of the parameters of trait-environment and migration-environment 
  # relationships, based on observed community composition 
  sink(tempfile())
  res.tab2 <- coalesc_abc(tab2, pool, multi = "tab", f.sumstats = f.sumstats,
                          filt.abc = function(x, par, env) filt_gaussian_env(x, par, env[1]), 
                          migr.abc = function(par, env) migr_env(par, env[2]),
                          par.filt = par.filt.range, par.migr = par.migr.range,
                          add = T, var.add = env,
                          nb.samp = nb.samp, parallel = FALSE, tol = 0.5,
                          pkg = NULL, #c("e1071","vegan"),
                          method.abc = "neuralnet")
  sink()
  expect_is(res.tab2, "list")
  expect_named(res.tab2, c("par", "obs", "obs.scaled", "ss", "ss.scale", "abc", "call"))
  expect_is(res.tab2$par, "matrix")
  expect_is(res.tab2$obs, "array")
  expect_is(res.tab2$obs.scaled, "array")
  expect_is(res.tab2$ss, "matrix")
  expect_is(res.tab2$ss.scale, "data.frame")
  expect_is(res.tab2$abc, "abc")
  sink()
  
  # 3/ Distinct pools of migrants across communities
  
  # The immigrants can be drawn from distinct pools for each community, to represent a local context.
  # In this case, the pool argument sent to coalesc_abc must be a list of the local pools.
  
  # Simulate several communities with distinct pools and same environmental filter
  pool.loc <- c()
  meta3 <- c()
  theta <- 2.5*(1:20) # The pools differ in diversity
  # We first define a baseline pool of immigrants common to all communities
  baseline <- coalesc(J = 10000, theta = 25)$com
  baseline$sp <- paste(0, baseline$sp, sep = "_")
  for(i in 1:20)
  {
    pool.loc[[i]] <- coalesc(J = 10000, theta = theta[i])$com
    pool.loc[[i]]$sp <- paste(i, pool.loc[[i]]$sp, sep = "_")
    pool.loc[[i]] <- rbind(pool.loc[[i]], baseline)
    meta3[[i]] <- coalesc(J = 400, pool = pool.loc[[i]], m = 0.5, 
                          filt = function(x) filt_gaussian(x, c(0.25, 0.1)))
  }
  # Build a species-by-site table
  tab3 <- array(0, c(20, length(unique(Reduce(rbind, pool.loc)[,2]))))
  colnames(tab3) <- unique(Reduce(rbind, pool.loc)[,2])
  for(i in 1:20) {
    compo <- table(meta3[[i]]$com$sp)
    tab3[i, names(compo)] <- compo
  }
  tab3 <- tab3[, colSums(tab3)!=0]
  
  # ABC estimation of the parameters of trait-environment and migration-environment 
  # relationships, based on observed community composition 
  ## Warning: this function may take a while
  par.range <- data.frame(rbind(c(0, 1), c(0.05, 1)))
  row.names(par.range) <- c("topt", "sigmaopt")
  sink(tempfile())
  res.tab3 <- coalesc_abc(tab3, pool.loc, multi = "tab", f.sumstats = f.sumstats,
                          filt.abc = filt_gaussian,
                          par.filt = par.range, 
                          nb.samp = nb.samp, parallel = FALSE, tol = 0.5,
                          pkg = NULL, #c("e1071","vegan"),
                          method.abc = "neuralnet")
  sink()
  expect_is(res.tab3, "list")
  expect_named(res.tab3, c("par", "obs", "obs.scaled", "ss", "ss.scale", "abc", "call"))
  expect_is(res.tab3$par, "matrix")
  expect_is(res.tab3$obs, "array")
  expect_is(res.tab3$obs.scaled, "array")
  expect_is(res.tab3$ss, "matrix")
  expect_is(res.tab3$ss.scale, "data.frame")
  expect_is(res.tab3$abc, "abc")
})