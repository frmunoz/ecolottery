# Returns a data frame of two variables which correlate with a population correlation of rho
# If desired, one of both variables can be fixed to an existing variable by specifying x

tcor <- function(n, rho = 0.5, mar.fun = rnorm, x = NULL, ...) {
  if (!is.null(x)) {X1 <- x} else {X1 <- mar.fun(n, ...)}
  if (!is.null(x) & length(x) != n) warning("Variable x does not have the same length as n!")
  
  if (rho > 1 | rho < 0) stop("rho must belong to [0; 1] value interval")
  
  C <- matrix(rho, nrow = 2, ncol = 2)
  diag(C) <- 1
  
  C <- chol(C)
  
  X2 <- mar.fun(n)
  X <- cbind(X1, X2)
  
  # induce correlation (does not change X1)
  df <- X %*% C
  
  ## if desired: check results
  #all.equal(X1, X[,1])
  #cor(X)
  
  # Formatting result into a dataframe
  df <- as.data.frame(df)
  # Naming columns
  if (is.null(x)){
    colnames(df) <- c("t1", "t2")
  } else{
    colnames(df) <- c("trait", "t2")
  }
  
  return(df)
}