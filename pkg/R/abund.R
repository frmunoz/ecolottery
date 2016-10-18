# From an object output by 'forward()' or 'coalesc()' computes vector of
# abundances per species
abund <- function(x)
{
  if(!is.null(x$com)) {
    # If there is a local pool computes regional level abundances
    loc <- as.data.frame(table(x$com[, "sp"]))
    colnames(loc) <- c("sp", "ab")
    loc$abrel <- loc$ab / length(x$com$ind)
  } else {
    stop("No local community provided, generate it with coalesc function")
  }
  
  if(!is.null(x$pool)) {
    reg <- as.data.frame(table(x$pool[, "sp"]))
    colnames(reg) <- c("sp", "ab")
    reg$abrel <- reg$ab / length(x$pool$ind)
    # If there is only a local community
    return(list(com = loc, pool = reg))
  } else {
    warning("No regional pool provided")
    return(list(com = loc))
  }
}
