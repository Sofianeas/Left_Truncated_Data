# Estimateur de régression de Nadaraya-Watson corrigé
NW_corrected <- function(x, X, Y, G_hat, h, K) {
  num <- sum((1 / G_hat) * Y * K((x - X) / h))
  den <- sum((1 / G_hat) * K((x - X) / h))
  if (den == 0) return(NA)
  return(num / den)
}
