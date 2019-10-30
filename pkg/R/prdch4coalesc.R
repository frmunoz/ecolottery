prdch4coalesc <- function(com.obs, pool, filt, params, stats = "abund", f.stats = NULL, estim = NULL){
  require(ggplot2)
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  sp.rank <- function(com){
    temp <- sort(table(com[,2]), decreasing = T)
    return(data.frame(sp=names(temp), rank=1:length(temp), abund = sapply(temp, function(x) x /sum(temp)), stringsAsFactors = F))
  }
  
  rel.abund <- function(com){
    temp <- sort(table(com[,2]), decreasing = T)
    sapply(temp, function(x) x /sum(temp))
  }
  
  if("custom" %in% stats) if(is.null(f.stats)) stop("Users should provide a function for computing summary statistics")
  if(ncol(data.frame(params))>1) if(is.null(filt)) stop("When simulating neutral communities params should only contain the migration rate posterior parameter distribution")
  if(ncol(data.frame(params))==1) if(!is.null(filt)) stop("Users should provided more than one posterior parameter distribution when simulating communities undergoing environmental filtering")
  
  lim <- 100 #à définir en fonction du temps necessaire à faire tourner ce nombre de simulation, possible de mettre une option pour forcer?
  J <- nrow(com.obs)
  
  if(!is.null(estim)){
    params <- switch(estim, 
                     "mean" = apply(params,2,mean),
                     "median" = apply(params,2,median),
                     "mode" = apply(params,2,getmode))
    
    params <- matrix(rep(params,each=lim),nrow=lim)
  }else{
    if(nrow(params) > lim){
      warning("Parameter distribution is large - Predictive checks using a sample of the posterior parameter distribution (ecrire ça mieux)")
      params <- params[sample(nrow(params),lim),]
    }
  }
  
  all.sims <- list()
  pb <- txtProgressBar(min = 0, max = nrow(params), style = 3)
  for(i in seq(nrow(params))){
    if(!is.null(filt)){
      par <- params[i,]
      
      sim <- ecolottery::coalesc(J, par[length(par)], theta = NULL, filt = function(x) filt(x, par[-length(par)]), 
                                 add = F,var.add =NULL, pool = pool, traits = NULL, 
                                 Jpool = 50 * J, 
                                 verbose = FALSE)}else{
      sim <-  ecolottery::coalesc(J, par[length(par)], theta = NULL, filt = NULL, 
                                  add = F,var.add =NULL, pool = pool, traits = NULL, 
                                  Jpool = 50 * J, 
                                  verbose = FALSE)}
                                  
                                  
                                
    all.sims[[i]] <- sim$com
    setTxtProgressBar(pb, i)
  }
  ret <- list()
  
  if("custom" %in% stats){
    ss.obs <- f.stats(com.obs)
    ss.sim <- lapply(all.sims, function(x) f.stats(x))
    ss.sim <- data.frame(matrix(unlist(ss.sim), nrow=length(ss.sim), byrow=T))
    pval <- sapply(seq(ncol(ss.sim)), function(x) sum(ss.sim[,x]<ss.obs[x])/(length(ss.sim[,x])+1))
    names(pval) <- paste0("S",seq(length(pval)))
    
    ret$pvalue <- pval
  }
  
  if("moments" %in% stats){
    moments <- function(com) moments::all.moments(as.numeric(com[,3]), order.max = 4)[-1]
    ss.obs <- moments(com.obs)
    ss.sim <- lapply(all.sims, function(x) moments(x))
    ss.sim <- data.frame(matrix(unlist(ss.sim), nrow=length(ss.sim), byrow=T))
    mval <- sapply(seq(ncol(ss.sim)), function(x) sum(ss.sim[,x]<ss.obs[x])/(length(ss.sim[,x])+1))
    names(mval) <- paste0("M",seq(length(mval)))
    
    ret$prd.moments <- mval
  }
  
  if("abund" %in% stats){
    #comparing relative abundances:
    abund.obs <- rel.abund(com.obs)
    abund.sim <- unlist(lapply(all.sims, function(x) rel.abund(x)))
    abund.sim <- split(abund.sim, names(abund.sim))
    spval <- sapply(seq(length(abund.obs)), FUN = function(i) sum(abund.sim[[names(abund.obs[i])]]<abund.obs[i])/(length(abund.sim[[names(abund.obs[i])]])+1))
    names(spval) <- names(abund.obs)
    
    
    
    #comparing rank-abundance curves:
    rank.obs <- cbind(sp.rank(com.obs), mdl=rep("OK", length(unique(com.obs$sp))),stringsAsFactors = FALSE)
    rank.obs[names(spval[which(spval<0.05)]),"mdl"] <- "Under-estimated"
    rank.obs[names(spval[which(spval>0.95)]),"mdl"] <- "Over-estimated"
    
    rank.sim <- dplyr::bind_rows(lapply(all.sims, function(x) sp.rank(x)))
    #dat <- rbind(cbind(rank.sim, type = rep("simulated", nrow(rank.sim))), cbind(rank.obs[,1:3], type = rep("observed", nrow(rank.obs))))
    
    
    d <- data.frame(rank = seq(length(unique((rank.sim$rank)))),
                    mean = unlist(lapply(split(rank.sim, as.factor(rank.sim$rank)), function(x) mean(x$abund))),
                    sd = unlist(lapply(split(rank.sim, as.factor(rank.sim$rank)), function(x) sd(x$abund))))
    
    cols <- c("Simulated"="#E69F00","Observed"="#56B4E9")
    print(ggplot(d, aes(x=rank, y=mean)) + 
      geom_point(data = rank.obs, aes(x=rank, y=abund, color="Observed")) +
      geom_line(data = rank.obs, aes(x=rank, y=abund, color="Observed"))+
      labs(x = "Rank", y = "Relative Abundance", title = "Species Rank-Abundance Curve")+
      scale_colour_manual(name=" ",values=cols) +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color ="Simulated"), width=0,position=position_dodge(0.05)))
    
    print(ggplot(rank.obs, aes(x=rank, y=abund, color=mdl, label=sp)) +
      geom_point(aes(color=mdl))+
      geom_text(angle = 90,hjust=0, vjust=0, size=3, show.legend = F)+
      labs(x = "Rank", y = "Relative Abundance", title = "Species Rank-Abundance Curve")+
      scale_color_manual(breaks = c("Under-estimated", "Over-estimated", "OK"),
                         values=c("#999999", "#E69F00", "#56B4E9")))
    
    ret$underestim.sp <- names(spval[which(spval<0.05)])
    ret$overestim.sp <- names(spval[which(spval>0.95)])
    
  }
  close(pb)
  return(ret)
}
