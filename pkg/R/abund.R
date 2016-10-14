abund <- function(x)
{
  if(!is.null(x$com)){
    loc <- as.data.frame(table(x$com[, "sp"]))
    colnames(loc) <- c("sp", "ab")
    loc$abrel <- loc$ab / length(x$com$ind)
  }
  else{
    stop("No local community provided, generate it with coalesc function")
  }
  
  if(!is.null(x$pool))
  {
    reg <- as.data.frame(table(x$pool[, "sp"]))
    colnames(reg) <- c("sp", "ab")
    reg$abrel <- reg$ab / length(x$pool$ind)
    return(list(loc = loc, reg = reg))
  } else{
    warning("No regional pool provided")
    return(list(loc = loc))
  }
}
