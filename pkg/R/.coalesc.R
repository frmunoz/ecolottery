# Simplified version of coalesc for faster computation in intensive calculations
.coalesc <- function (J, m = 1) 
{
  if (m < 1)
  {
    I <- m * (J - 1)/(1 - m)
    X <- runif(J)
    R_com <- I/(I + (1:J) - 1)
    assign_com <- which(X <= R_com)
    unassign_com <- which(X > R_com)
  } else {
    assign_com = 1:J
    unassign_com = NULL
  }
  com_species <- length(assign_com)
  migrants <- sample(1:nrow(pool), com_species, prob = filt(data.frame(pool[,-(1:2)])))
  ind_com_lab <- array(NA, J)
  ind_com_lab[assign_com] <- pool[migrants[1:com_species], 1] # Assign individuals
  ind_com_sp <- array(NA, J)
  ind_com_sp[assign_com] <- pool[migrants[1:com_species], 2] # Assign species
  ind_com_traits <- matrix(NA,nrow=J,ncol=ncol(data.frame(pool[,-(1:2)])))
  ind_com_traits[assign_com,] <- apply(data.frame(pool[,-(1:2)]),2,function(y) y[migrants[1:com_species]]) # Assign traits
  if (!is.null(unassign_com)) {
    for (j in 1:length(unassign_com)) {
        existing_sp <- sample.int(unassign_com[j] - 1, 1)
        ind_com_lab[unassign_com[j]] <- ind_com_lab[existing_sp]
        ind_com_sp[unassign_com[j]] <- ind_com_sp[existing_sp]
        ind_com_traits[unassign_com[j], ] <- ind_com_traits[existing_sp, ]
      }
  }
  com <- data.frame(ind = ind_com_lab, sp = ind_com_sp, ind_com_traits)
  if (m == 1 & is.null(filt)) return(list(pool = com)) else return(list(com = com, pool = pool))
}
