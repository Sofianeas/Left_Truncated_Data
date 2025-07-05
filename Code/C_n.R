# Estimation de la fonction de troncature C_n(y)
C_n <- function(y, T, Y) {
  mean(T <= y & y <= Y)
}