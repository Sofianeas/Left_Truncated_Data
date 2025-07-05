# Fonction pour entraîner un réseau de neurones avec correction de troncature
train_nn_estimator <- function(X_obs, Y_obs, T_obs, 
                               size = 10, decay = 0.01, 
                               maxit = 1000, linout = TRUE,
                               x_grid = seq(0, 1, length.out = 100)) {
  
  # 1. Calculer les poids de troncature via Lynden-Bell
  n <- length(X_obs)
  GG <- est_GLBB(Y_obs, T_obs, n)
  weights <- 1 / pmax(GG, 1e-6)  # Éviter la division par zéro
  
  # 2. Préparer les données
  train_data <- data.frame(X = X_obs, Y = Y_obs)
  
  # 3. Entraîner le réseau de neurones avec les poids
  set.seed(123)  # Pour la reproductibilité
  nn_model <- nnet::nnet(Y ~ X, 
                         data = train_data, 
                         weights = weights,
                         size = size,        # Nombre d'unités cachées
                         decay = decay,      # Paramètre de régularisation
                         maxit = maxit,      # Itérations maximales
                         linout = linout,    # Sortie linéaire pour la régression
                         trace = FALSE)      # Désactive les logs
  
  # 4. Faire des prédictions sur la grille
  predictions <- predict(nn_model, newdata = data.frame(X = x_grid))
  
  list(
    model = nn_model,
    predictions = predictions,
    x_grid = x_grid
  )
}

# Fonction pour comparer avec l'estimateur NW
compare_nn_vs_nw <- function(n, m_true, T_dist, T_params, 
                             nn_params, h = 0.1, n_sim = 10, seed = NULL) {
  
  mse_nn <- numeric(n_sim)
  mse_nw <- numeric(n_sim)
  x_grid <- seq(0.05, 0.95, length.out = 100)
  true_vals <- m_true(x_grid)
  
  for (i in 1:n_sim) {
    if (!is.null(seed)) set.seed(seed + i)
    data <- simulate_data(n, m_true, T_dist, T_params)
    
    # Estimateur NW
    nw_est <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h)
    mse_nw[i] <- mean((nw_est - true_vals)^2, na.rm = TRUE)
    
    # Réseau de neurones
    nn_result <- train_nn_estimator(
      X_obs = data$X_obs,
      Y_obs = data$Y_obs,
      T_obs = data$T_obs,
      size = nn_params$size,
      decay = nn_params$decay,
      maxit = nn_params$maxit,
      x_grid = x_grid
    )
    mse_nn[i] <- mean((nn_result$predictions - true_vals)^2, na.rm = TRUE)
  }
  
  list(
    mse_nw = mean(mse_nw),
    mse_nn = mean(mse_nn),
    efficiency_ratio = mean(mse_nn) / mean(mse_nw)
  )
}