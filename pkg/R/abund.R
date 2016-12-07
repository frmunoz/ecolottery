# From an object output by 'forward()' or 'coalesc()' computes vector of
# abundances per species
abund <- function(x) {
  if (!is.list(x)) {
    stop("Provided argument needs to be a list of communities")
  }
  
    rel_abund_list <-  lapply(x, function(y) {
      
      if (is.list(y)) {
        rel_abund <- lapply(y, get_rel_abund)
      } else {
        rel_abund <- get_rel_abund(y)
      } 
      
      return(rel_abund)
    })
    
    if (sum(is.na(rel_abund_list)) + sum(is.null(rel_abund_list)) != 0) {
      warning("Some communities were undefined; returning NA abundances")
    }
    
    return(rel_abund_list)
}

# Internal function to compute relative abundances from a community data.frame
get_rel_abund <- function(comdf) {
  
  if (is.data.frame(comdf) | is.matrix(comdf)) {
    rel_abund <- as.data.frame(table(comdf[, "sp"]))
    colnames(rel_abund) <- c("sp", "ab")
    rel_abund$relab <- rel_abund$ab / nrow(comdf)
  } else {
    rel_abund <- NA
  }
  
  return(rel_abund)
}