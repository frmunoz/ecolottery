# From an output object of 'forward()' or 'coalesc()', computes vector of
# abundances per species
abund <- function(x) {
  if (!is.list(x)) {
    stop("The input argument must be a list of assemblages", call. = FALSE)
  }
  
  # Select only elements that are data frames
  if("com_t" %in% names(x)) {
    y <- x$com_t[lapply(x$com_t, class)=="data.frame"]
    names(y) <- paste("com_t", 1:length(x$com_t), sep=".")
  } else y <- NULL
  x <- x[lapply(x, class)=="data.frame"]
  x <- c(x, y)
  
  rel_abund_list <-  lapply(x, function(y) {
      
      if (is.list(y) & !is.data.frame(y) & !is.matrix(y)) {
        rel_abund <- lapply(y, .get_rel_abund)
      } else {
        rel_abund <- .get_rel_abund(y)
      } 
      
      return(rel_abund)
  })
    
  if (sum(is.na(rel_abund_list)) + sum(is.null(rel_abund_list)) != 0) {
      warning("Some communities were undefined; returning NA abundances",
              call. = FALSE)
  }
  
  return(rel_abund_list)
}

# Internal function to compute relative abundances from a community data.frame
.get_rel_abund <- function(comdf) {
  
  if (is.data.frame(comdf) | is.matrix(comdf)) {
    rel_abund <- as.data.frame(table(comdf[, "sp"]), stringsAsFactors = FALSE)
    rownames(rel_abund) <- rel_abund[,1]
    rel_abund$relab <- rel_abund[,2] / nrow(comdf)
    rel_abund <- rel_abund[,-1]
    colnames(rel_abund)[1] <- "ab"
  } else {
    rel_abund <- NA
  }
  
  return(rel_abund)
}
