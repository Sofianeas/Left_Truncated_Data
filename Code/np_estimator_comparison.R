comparison_results <- eventReactive(input$run_comparison, {
  req(input$comp_n, input$comp_h)
  
  m_true <- parse_true_func()
  T_params <- get_t_params()
  seed <- if (is.na(input$seed)) NULL else input$seed
  
  withProgress(message = 'Running estimator comparison', value = 0, {
    # Simulate data
    set.seed(seed)
    data <- simulate_data(input$comp_n, m_true, input$t_dist, T_params)
    x_grid <- seq(0.05, 0.95, length.out = 100)
    m_vals <- m_true(x_grid)
    
    # Prepare results list
    results <- list()
    
    # 1. Nadaraya-Watson (corrected for truncation)
    incProgress(0.15, detail = "Running NW estimator")
    nw_est <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, input$comp_h)
    results$nw <- list(estimate = nw_est, mse = mean((nw_est - m_vals)^2, na.rm = TRUE))
    
    # 2. LOESS
    incProgress(0.15, detail = "Running LOESS")
    df_loess_data <- data.frame(x = data$X_obs, y = data$Y_obs)
    loess_est <- tryCatch({
      predict(loess(y ~ x, data = df_loess_data, span = input$comp_loess_span), 
              newdata = data.frame(x = x_grid))
    }, error = function(e) rep(NA, length(x_grid)))
    results$loess <- list(estimate = loess_est, mse = mean((loess_est - m_vals)^2, na.rm = TRUE))
    
    # 3. Local polynomial
    incProgress(0.15, detail = "Running local polynomial")
    local_poly_est <- tryCatch({
      fit_poly <- locfit(y ~ lp(x, deg = input$comp_poly_deg, h = input$comp_h), 
                         data = df_loess_data)
      predict(fit_poly, newdata = x_grid)
    }, error = function(e) rep(NA, length(x_grid)))
    results$poly <- list(estimate = local_poly_est, mse = mean((local_poly_est - m_vals)^2, na.rm = TRUE))
    
    # 4. Splines
    incProgress(0.15, detail = "Running splines")
    spline_est <- tryCatch({
      predict(smooth.spline(data$X_obs, data$Y_obs, spar = input$comp_spline_spar), 
              x = x_grid)$y
    }, error = function(e) rep(NA, length(x_grid)))
    results$spline <- list(estimate = spline_est, mse = mean((spline_est - m_vals)^2, na.rm = TRUE))
    
    # 5. Neural Network (nouvelle partie à ajouter)
    incProgress(0.2, detail = "Running neural network")
    nn_est <- tryCatch({
      nn_result <- train_nn_estimator(
        X_obs = data$X_obs,
        Y_obs = data$Y_obs,
        T_obs = data$T_obs,
        layers = input$nn_layers,
        units = input$nn_units,
        epochs = input$nn_epochs,
        batch_size = input$nn_batch,
        activation = input$nn_activation,
        lr = input$nn_lr,
        x_grid = x_grid
      )
      nn_result$predictions
    }, error = function(e) rep(NA, length(x_grid)))
    results$nn <- list(estimate = nn_est, mse = mean((nn_est - m_vals)^2, na.rm = TRUE))
    
    incProgress(0.2, detail = "Finalizing results")
    
    list(
      x_grid = x_grid,
      m_true = m_vals,
      results = results,
      data = data
    )
  })
})



