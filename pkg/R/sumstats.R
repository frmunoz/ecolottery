sumstats <- function(com, multi = "single", traits = NULL, type = "mix", n = 4, m = 4)
{
  if (!requireNamespace("hillR", quietly = TRUE)) stop("sumstats requires package hillR to be installed")
  if (!requireNamespace("moments", quietly = TRUE)) stop("function requires package moments to be installed")
  
  if(!type%in%c("taxo","func","mix")) stop("type should be either taxo, func or mix")
  if(!multi%in%c("single","tab")) stop("mutli should be single or tab")
  
  
  if(multi == "single"){
    #if(any(is.na(com[,c(1,2)]))) stop("NA in community")
    #if(ncol(com) > 3) stop("community data.frame/matrix should have at most 3 columns: individuals, species and trait values (in that order)")
    
    if(type == "taxo"){
      if((ncol(com) < 2) & (!is.null(traits))) warning("trait information is ignored when type is taxo")
      tab <- table(com[,2])
      Ss <- as.vector(t(sapply(0:(n-1), function(x) hillR::hill_taxa(tab, q=x))))
    }
    
    if(type == "func"){
      if(ncol(com) < 3){
        if(is.null(traits)){
          stop("trait information must be provided")
        }else{
          if(any(is.na(traits))){ warning("Individuals with missing trait information have been removed")
            com <- com[-which(com[,2] == rownames(traits)[which(is.na(traits))]),] }
          
          com[,3] <- traits[as.character(com[,2]),]} } 
      if(any(is.na(com[,3]))) { warning("Individuals with missing trait information have been removed"); com <- com[complete.cases(com),]}
      Ss <- moments::all.moments(com[,3], order.max = m)[-1]
      
    }
    if(type == "mix"){
      if(ncol(com)<3){
        if(is.null(traits)){
          stop("trait information must be provided")
        }else{if(any(is.na(traits))) { warning("Individuals with missing trait information have been removed"); com <- com[-which(com[,2] == rownames(traits)[which(is.na(traits))]),] }
          com[,3] <- traits[as.character(com[,2]),]} 
      }else{ if(any(is.na(com[,3]))) {warning("Individuals with missing trait information have been removed")
        com <- com[complete.cases(com),] } } 
      
      tab <- table(com[,2])
      
      Ss <- c(moments::all.moments(as.numeric(com[,3]), order.max = m)[-1],
              as.vector(t(sapply(0:(n-1), function(x) hillR::hill_taxa(tab, q=x)))))
      
    }
  }
  
  if(multi == "tab"){
    if(any(is.na(com))) stop("NA in species by site matrix/data.frame")
    if(type == "taxo"){
      Ss <- as.vector(t(sapply(0:(n-1), function(q) hillR::hill_taxa(com,q))))
      
    }
    if(type == "func"){
      if(is.null(traits))stop("trait information should be provided in traits")
      if(any(is.na(traits))){ warning("Species with missing trait information have been removed from species by site matrix/data.frame") 
        com <- com[,!colnames(com)==rownames(traits)[which(is.na(traits))]] }
      
      Ss <- as.vector(apply(com, 1, FUN = function(x) moments::all.moments(traits[rep(colnames(com),x),], order.max = m)[-1]))
    }
    
    if(type == "mix"){
      if(is.null(traits)) stop("trait information should be provided in traits")
      if(any(is.na(traits))){ warning("Species with missing trait information have been removed from species by site matrix/data.frame") 
        com <- com[,!colnames(com)==rownames(traits)[which(is.na(traits))]] }
      
      Ss <- c(as.vector(apply(com, 1, FUN = function(x) moments::all.moments(traits[rep(colnames(com),x),][!is.na(traits[rep(colnames(com),x),])], order.max = m)[-1])),
              as.vector(t(sapply(0:(n-1), function(q) hillR::hill_taxa(com,q)))))
      
    }
  }
  
  return(Ss)
}
