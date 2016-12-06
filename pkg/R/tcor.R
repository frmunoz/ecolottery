# Generates a data.frame with two correlated traits with the given 'mar.fun'
# function correlated by a factor 'rho'. Or generates a trait vector correlated
# to the given one 'x' by a factor 'rho'. 'n' gives the number of values

tcor <- function(n, rho = 0.5, mar.fun = rnorm, x = NULL, ...) {
  
  
  if (!is.null(x)) {
    # If the user provides a target trait vector uses it
    first_trait <- x
  } else {
    # Otherwise generates the first trait vector with provided function
    first_trait <- mar.fun(n, ...)
  }
  
  # If the provided vector by the user is different from the target vector size
  if (!is.null(x) & length(x) != n) {
    warning("Trait vector x does not have length n!")
  }
  
  # correlation needs to be strictly between 0 and 1 
  if (rho > 1 | rho < 0) {
    stop("rho must belong to ]0; 1[ value interval")
  }
  
  # Generate a correlation matrix with given correlation coefficient
  corr_mat <- matrix(rho, nrow = 2, ncol = 2)
  diag(corr_mat) <- 1
  
  corr_mat <- chol(corr_mat)
  
  second_trait <- mar.fun(n)
  trait_mat <- cbind(first_trait, second_trait)
  
  # Induces correlation (does not change first_trait)
  trait_df <- trait_mat %*% corr_mat
  
  # Formatting result into a dataframe
  trait_df <- as.data.frame(trait_df)
  # Naming columns
  if (is.null(x)) {
    colnames(trait_df) <- c("t1", "t2")
  } else{
    colnames(trait_df) <- c("trait", "t2")
  }
  
  return(trait_df)
}