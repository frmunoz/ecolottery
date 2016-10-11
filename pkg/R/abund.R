# From an object output by 'forward()' or 'coalesc()' computes vector of
# abundances per species
abund <- function(x)
{
  # Abundances of species in local community
  loc_t <- as.data.frame(table(x$com[,2]))
  loc <- loc_t[,2]  # Get a flat vector of abundances
  names(loc) <- loc_t[,1]
  
  if (!is.null(x$pool))
  {
    # If the is a local pool computes regional level abundances
    reg_t <- as.data.frame(table(x$pool[,2]))
    reg <- reg_t[,2]
    names(reg) <- reg_t[,1]
    return(list(com = loc, pool = reg))
  } else {
    # If there is only a local community
    return(list(com = loc))
  }
}
