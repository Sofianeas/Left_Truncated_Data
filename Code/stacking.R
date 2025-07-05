stacking_results <- eventReactive(input$run_stack_plot, {
  req(input$stack_n, input$stack_methods, input$stacking_method)
  
  m_true <- parse_true_func()
  T_params <- get_t_params()
  x_grid <- seq(0, 1, length.out = 100)
  h_grid <- seq(0.03, 0.4, length.out = 10)
  
  withProgress(message = "Simulation Stacking : courbes", value = 0, {
    
    data <- simulate_data(input$stack_n, m_true, input$t_dist, T_params)
    X_obs <- data$X_obs; Y_obs <- data$Y_obs; T_obs <- data$T_obs
    true_vals <- m_true(x_grid)
    
    # Résultats des estimateurs
    results <- list()
    df_meta <- data.frame(Y = Y_obs)
    df_pred <- data.frame(x = x_grid)
    
    # --- NW
    if ("nw" %in% input$stack_methods) {
      h_val <- if (isTRUE(input$auto_h_stack)) {
        select_h_cv(X_obs, Y_obs, T_obs, h_grid)
      } else input$stack_h
      m_nw <- nw_estimator_truncated(x_grid, X_obs, Y_obs, T_obs, h_val)
      results$nw <- m_nw
      df_meta$m_nw <- nw_estimator_truncated(X_obs, X_obs, Y_obs, T_obs, h_val)
      df_pred$m_nw <- m_nw
    }
    
    # --- NN
    if ("nn" %in% input$stack_methods) {
      nn_result <- train_nn_estimator(X_obs, Y_obs, T_obs,
                                      size = input$stack_nn_size,
                                      decay = input$stack_nn_decay,
                                      maxit = input$stack_nn_maxit,
                                      x_grid = x_grid)
      results$nn <- nn_result$predictions
      df_meta$m_nn <- predict(nn_result$model, newdata = data.frame(X = X_obs))
      df_pred$m_nn <- nn_result$predictions
    }
    
    # --- Apprentissage meta-learner choisi
    meta_model <- switch(input$stacking_method,
                         "lm" = lm(Y ~ ., data = df_meta),
                         "mars" = earth(Y ~ ., data = df_meta),
                         "svr" = svm(Y ~ ., data = df_meta),
                         "xgb" = {
                           dtrain <- xgb.DMatrix(data = as.matrix(df_meta[, -1]), 
                                                 label = df_meta$Y)
                           xgb.train(data = dtrain, nrounds = 100, 
                                     objective = "reg:squarederror", verbose = 0)
                         },
                         "rf" = randomForest(Y ~ ., data = df_meta),
                         NULL)
    
    # --- Prédiction stacking
    m_stack <- tryCatch({
      if (input$stacking_method == "xgb") {
        predict(meta_model, newdata = as.matrix(df_pred[, -1]))
      } else {
        predict(meta_model, newdata = df_pred[, -1, drop = FALSE])
      }
    }, error = function(e) rep(NA, length(x_grid)))
    
    incProgress(1)
    
    list(
      x = x_grid,
      m_true = true_vals,
      results = results,
      m_stack = m_stack,
      h_val = if(exists("h_val")) h_val else NULL,
      methods_used = input$stack_methods,
      meta_method = input$stacking_method
    )
  })
})