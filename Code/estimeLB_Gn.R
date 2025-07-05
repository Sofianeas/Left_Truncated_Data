# Estimateur de la fonction de répartition G à l'aide de Lynden-Bell

estimeLB_Gn <- function(T_obs, Y_obs) {
  n <- length(T_obs)
  G_LB <- numeric(n)
  
  for (i in seq_along(T_obs)) {
    indices <- which(T_obs > T_obs[i])
    product_terms <- (n * sapply(indices, function(j) mean(T_obs <= T_obs[j] & T_obs[j] <= Y_obs))) /
      (n * sapply(indices, function(j) mean(T_obs <= T_obs[j] & T_obs[j] <= Y_obs)) + 1)
    G_LB[i] <- if (length(product_terms) == 0) 1 else prod(product_terms)
  }
  
  return(G_LB)
}
