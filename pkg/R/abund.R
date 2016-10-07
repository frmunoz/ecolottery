abund <- function(x)
{
  loc_t <- as.data.frame(table(x$com[,2]))
  loc <- loc_t[,2]; names(loc) <- loc_t[,1]
  if(!is.null(x$pool))
  {
    reg_t <- as.data.frame(table(x$pool[,2]))
    reg <- reg_t[,2]; names(reg) <- reg_t[,1]
    return(list(loc=loc, reg=reg))
  } else return(list(loc=loc))
}
