# Estimateur de la densité de X
f_X_hat <- function(x, X, G_hat, h, K) {
  n <- length(X)
  alpha <- 1  # À adapter selon vos données
  (alpha / (n * h)) * sum((1 / G_hat) * K((x - X) / h))
}

# Estimateur de F1(x, y)
F1_hat <- function(x, y, X, Y, G_hat, hK, hH, K, Phi) {
  n <- length(X)
  alpha <- 1  # À adapter selon vos données
  (alpha / (n * hK)) * sum((1 / G_hat) * Phi((y - Y) / hH) * K((x - X) / hK))
}
