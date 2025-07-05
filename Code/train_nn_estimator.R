# Fonction pour entraîner un réseau de neurones avec correction de troncature
train_nn_estimator <- function(X_obs, Y_obs, T_obs, 
                               size = 10, decay = 0.01, 
                               maxit = 1000, linout = TRUE,
                               x_grid = seq(0, 1, length.out = 100)) {
  
  # Calculate truncation weights
  n <- length(X_obs)
  GG <- est_GLBB(Y_obs, T_obs, n)
  weights <- 1 / pmax(GG, 1e-6)  # Avoid division by zero
  
  # Prepare data
  train_data <- data.frame(X = X_obs, Y = Y_obs)
  
  # Train neural network
  set.seed(123)
  nn_model <- nnet(Y ~ X, data = train_data, weights = weights,
                   size = size, decay = decay, maxit = maxit,
                   linout = linout, trace = FALSE)
  
  # Make predictions
  predictions <- predict(nn_model, newdata = data.frame(X = x_grid))
  
  list(
    model = nn_model,
    predictions = predictions,
    x_grid = x_grid
  )
}

