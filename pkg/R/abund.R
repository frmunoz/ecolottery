# From an object output by 'forward()' or 'coalesc()' computes vector of
# abundances per species
abund <- function(x)
{
  if(!is.null(x$com)) {
    # If there is a local pool computes regional level abundances
    loc <- as.data.frame(table(x$com[, "sp"]))
    colnames(loc) <- c("sp", "ab")
    rownames(loc) <- loc$sp
    loc$relab <- loc$ab / nrow(x$com)
  } else {
    stop("No local community provided, generate it with coalesc function")
  }
  
  if(!is.null(x$pool)) {
    reg <- as.data.frame(table(x$pool[, "sp"]))
    colnames(reg) <- c("sp", "ab")
    rownames(reg) <- reg$sp
    reg$relab <- reg$ab / nrow(x$pool)
    # If there is only a local community
    return(list(com = loc, pool = reg))
  } else {
    warning("No regional pool provided")
    return(list(com = loc))
  }
}
