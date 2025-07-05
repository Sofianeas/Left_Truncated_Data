# Fonction de simulation de données tronquées
simulate_truncated_data <- function(n, m_function, error_sd = 1, truncation_dist = function(n) runif(n, min = 0, max = 2)) {
  N <- n * 2  # Génère plus de données pour compenser la troncature
  X <- runif(N, min = 0, max = 2)
  epsilon <- rnorm(N, mean = 0, sd = error_sd)
  Y <- m_function(X) + epsilon
  T <- truncation_dist(N)
  
  keep <- which(Y >= T)
  if (length(keep) < n) stop("Pas assez de données après troncature")
  indices <- sample(keep, n)
  
  list(X = X[indices], Y = Y[indices], T = T[indices])
}
