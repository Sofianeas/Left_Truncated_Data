# Brouillon App finale à présenter

library(shiny)
library(plotly)
library(shinythemes)
library(DT)
library(bslib)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(locfit)
library(markdown)
library(xfun)
library(nnet)
library(future.apply)
library(compiler)
library(compiler)
library(xgboost)
library(e1071)
library(randomForest)
library(earth)
library(tidyr)
library(dplyr)
plan(multisession)  # ou multicore si tu es sur Linux/macOS
enableJIT(3)

# Custom theme
my_theme <- bs_theme(
  version = 5,
  bg = "#f8f9fa",
  fg = "#212529",
  primary = "#2c3e50",
  secondary = "#6c757d",
  success = "#28a745",
  base_font = font_google("Roboto"),
  heading_font = font_google("Montserrat"),
  "font-size-base" = "1.1rem"
)

# Fonction C_n pour estimer P(T <= y <= Y)
C_n <- function(Y, T, y, n) {
  mean((T <= y) & (Y >= y))
}

# Estimateur de Lynden-Bell pour la fonction de survie de T
estimeLB_Gn <- function(Y, T, t, n) {
  prod <- 1
  for (i in 1:n) {
    if (T[i] > t) {
      c_val <- C_n(Y, T, T[i], n)
      prod <- prod * ((n * c_val) / (n * c_val + 1))
    }
  }
  return(prod)
}

# Version vectorisée de l'estimateur de Lynden-Bell
est_GLBB <- function(Y, T, n) {
  sapply(Y, function(y) estimeLB_Gn(Y, T, y, n))
}

# Noyau gaussien
noy <- function(x) {
  (1/sqrt(2*pi)) * exp((-1/2) * x^2)
}

# Fonction H (utilisant la normale standard)
fonction_H <- function(y) {
  pnorm(y, mean = 0, sd = 1)
}

# Estimateur de la densité de X
est_densityv <- function(X, hK, x, n, alpha, GG) {
  sum_val <- sum((1/GG) * noy((x - X)/hK))
  (alpha/(n * hK)) * sum_val
}

# Estimateur de F1(x,y)
est2F1 <- function(X, Y, hK, hH, x, y, n, alpha, GG) {
  sum_val <- sum(fonction_H((y - Y)/hH) * (1/GG) * noy((x - X)/hK))
  (alpha/(n * hK)) * sum_val
}

# Fonctions supplémentaires pour l'estimation de m(x) = E[Y|X=x]

# Estimateur NW corrigé pour données tronquées
nw_estimator_truncated <- function(x_grid, X_obs, Y_obs, T_obs, h) {
  n <- length(X_obs)
  
  # 1. Estimation Lynden-Bell une seule fois
  GG <- est_GLBB(Y_obs, T_obs, n)
  GG <- pmax(GG, 1e-6)
  
  # 2. Matrice des noyaux
  K_matrix <- outer(x_grid, X_obs, function(x, xi) dnorm((x - xi) / h))
  
  # 3. Pondération par 1 / GG_j (colonne par colonne)
  weights_mat <- sweep(K_matrix, 2, GG, FUN = "/")
  
  # 4. Estimation vectorisée
  numerator <- weights_mat %*% Y_obs
  denominator <- rowSums(weights_mat)
  
  # 5. Sécurité : éviter division par 0
  m_est <- ifelse(denominator < 1e-6, NA, numerator / denominator)
  
  return(as.vector(m_est))
}



simulate_data <- function(n, m_true, t_dist = "normal", t_params = c(0, 0.5), 
                          trunc_method = "random", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- runif(n, 0, 1)
  epsilon <- rnorm(n, sd = 0.2)
  Y <- m_true(X) + epsilon
  
  if (trunc_method == "fixed") {
    # Taux de troncature fixé
    target_rate <- t_params[1]
    Y_sorted <- sort(Y)
    threshold_idx <- max(1, floor(n * target_rate))
    T <- rep(Y_sorted[threshold_idx], n)
  } else {
    # Troncature aléatoire
    T <- switch(t_dist,
                "normal" = rnorm(n, mean = t_params[1], sd = t_params[2]),
                "uniform" = runif(n, min = t_params[1], max = t_params[2]),
                "exponential" = rexp(n, rate = t_params[1]),
                "gamma" = rgamma(n, shape = t_params[1], rate = t_params[2]),
                rnorm(n, mean = 0, sd = 0.5)) # défaut
  }
  
  idx <- which(Y >= T)
  if (length(idx) == 0) idx <- which.max(Y - T)
  
  list(
    X_obs = X[idx], 
    Y_obs = Y[idx], 
    T_obs = T[idx],
    X_full = X,
    Y_full = Y,
    T_full = T
  )
}



select_h_cv <- function(X, Y, T, h_grid, K = 3) {
  n <- length(Y)
  folds <- sample(rep(1:K, length.out = n))
  GG <- est_GLBB(Y, T, n)
  GG <- pmax(GG, 1e-6)
  
  cv_errors <- sapply(h_grid, function(h) {
    errors <- numeric(K)
    for (k in 1:K) {
      idx_train <- which(folds != k)
      idx_test <- which(folds == k)
      x_train <- X[idx_train]
      y_train <- Y[idx_train]
      t_train <- T[idx_train]
      x_test <- X[idx_test]
      y_test <- Y[idx_test]
      g_test <- GG[idx_test]
      m_hat <- nw_estimator_truncated(x_test, x_train, y_train, t_train, h)
      errors[k] <- mean(((y_test - m_hat)^2) / g_test, na.rm = TRUE)
    }
    mean(errors)
  })
  
  h_grid[which.min(cv_errors)]
}





# Analysis function for bias, variance, MSE
analyze_performance <- function(n, N_sim, m_true, T_dist, T_params,
                                h_grid = seq(0.03, 0.4, length.out = 10),
                                h_fixed = NULL, seed = NULL) {
  x_grid <- seq(0.1, 0.9, length.out = 20)
  m_hat_mat <- matrix(NA, nrow = N_sim, ncol = length(x_grid))
  h_used <- numeric(N_sim)
  
  for (b in 1:N_sim) {
    if (!is.null(seed)) set.seed(seed + b)
    
    data <- simulate_data(n, m_true, T_dist, T_params)
    
    h_opt <- if (is.null(h_fixed)) {
      select_h_cv(data$X_obs, data$Y_obs, data$T_obs, h_grid)
    } else {
      h_fixed
    }
    
    h_used[b] <- h_opt
    
    m_hat <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h_opt)
    m_hat_mat[b, ] <- m_hat
  }
  
  m_mean <- colMeans(m_hat_mat, na.rm = TRUE)
  m_var <- apply(m_hat_mat, 2, var, na.rm = TRUE)
  bias_sq <- (m_mean - m_true(x_grid))^2
  mse <- bias_sq + m_var
  
  list(
    x = x_grid,
    bias = sqrt(bias_sq),
    variance = sqrt(m_var),
    mse = sqrt(mse),
    n = n,
    h_mean = mean(h_used),
    h_all = h_used,
    N_sim = N_sim
  )
}






# Function to compute MISE
compute_mise <- function(N, m_true, T_dist, T_params,
                         N_sim = 50, grid_points = 50,
                         h_grid = seq(0.03, 0.4, length.out = 10),
                         seed = NULL) {
  
  grid_x <- seq(0, 1, length.out = grid_points)
  m_hat_mat <- matrix(0, nrow = N_sim, ncol = grid_points)
  h_used <- numeric(N_sim)
  
  for (i in 1:N_sim) {
    if (!is.null(seed)) set.seed(seed + i)
    
    data <- simulate_data(N, m_true, T_dist, T_params)
    h_opt <- select_h_cv(data$X_obs, data$Y_obs, data$T_obs, h_grid)
    h_used[i] <- h_opt
    
    m_hat <- nw_estimator_truncated(grid_x, data$X_obs, data$Y_obs, data$T_obs, h_opt)
    m_hat_mat[i, ] <- m_hat
  }
  
  true_vals <- m_true(grid_x)
  mise <- mean(rowMeans((m_hat_mat - matrix(true_vals, nrow = N_sim, ncol = grid_points, byrow = TRUE))^2, na.rm = TRUE))
  
  list(
    n = N,
    mise = mise,
    h_mean = mean(h_used),
    h_all = h_used
  )
}






# UI
ui <- page_navbar(
  title = "Nonparametric Estimation with Left-Truncated Data",
  theme = my_theme,
  nav_panel(
    "Theory & Code",
    card(
      card_header("Méthodologie & Code R", class = "bg-primary text-white"),
      card_body(
        withMathJax(HTML("
        <h4><strong>📘 Partie Théorique</strong></h4>

        

          <h5>1. Estimateur de Lynden-Bell \\( \\widehat{G}_{LB}(t) \\)</h5>
          <p>
            Cet estimateur adapte l'estimateur de Kaplan-Meier pour les données tronquées à gauche :
            <br>
            \\[ \\widehat{G}_{LB}(t) = \\prod_{i: T_i > t} \\left( \\frac{n \\widehat{C}_n(T_i)}{n \\widehat{C}_n(T_i) + 1} \\right) \\]
            Il estime la fonction de répartition \\(G(t) = \\mathbb{P}(T \\leq t)\\).
          </p>
          
          <h5>2. Estimation de la fonction C(y)</h5>
          <p>
            L'estimateur empirique de la fonction C(y) défini par :
            <br>
            \\[ \\widehat{C}_n(y) = \\frac{1}{n} \\sum_{i=1}^n \\mathbf{1}_{\\{T_i \\leq y \\leq Y_i\\}} \\]
            Il mesure la probabilité que l'observation soit visible au niveau \\(y\\).
          </p>

          <h5>3. Estimateur de la densité de \\(X\\)</h5>
          <p>
            L’estimateur de la densité marginale de \\(X\\) à l’aide du noyau gaussien est donné par :
            <br>
            \\[ \\widehat{f}_X(x) = \\frac{\\alpha}{n h_K} \\sum_{j=1}^n \\frac{1}{\\widehat{G}_j} K\\left(\\frac{x - X_j}{h_K}\\right) \\]
            où \\(K\\) est le noyau gaussien et \\(\\widehat{G} = \\widehat{G}_{LB}(Y_j)\\).
          </p>


          <h5>4. Estimateur de la fonction de régression \\( m(x) = \\mathbb{E}[Y | X = x] \\)</h5>
          <p>
            La forme pratique de l’estimateur de Nadaraya-Watson corrigé pour troncature est :
            <br>
            \\[
            \\widehat{m}_n(x) = 
            \\frac{ \\sum_{j=1}^n \\frac{1}{\\widehat{G}_{Y_{j}}} Y_j K\\left(\\frac{x - X_j}{h_K} \\right) }
                 { \\sum_{j=1}^n \\frac{1}{\\widehat{G}_{Y_{j}}} K\\left(\\frac{x - X_j}{h_K} \\right) }
            \\]
            Cette formule ajuste les poids des observations en fonction de la troncature.
          </p>

          <h5>5. Simulation de données tronquées</h5>
          <p>
            Pour évaluer les performances de l'estimateur \\( \\widehat{m}_n(x) \\), on simule des données selon le modèle :
            <br><br>
            \\[
            Y_i = m(X_i) + \\varepsilon_i, \\quad \\text{où} \\quad \\varepsilon_i \\sim \\mathcal{N}(0, \\sigma^2)
            \\]
            <br>
            Ensuite, on génère la variable de troncature \\(T_i\\) suivant une distribution choisie (normale, uniforme, exponentielle, gamma).
            L'observation \\((X_i, Y_i, T_i)\\) est conservée uniquement si \\(Y_i \\geq T_i\\), ce qui introduit une troncature à gauche.
          </p>
          <p>
            La fonction R <code>simulate_data()</code> permet de :
            <ul>
              <li>générer les variables \\(X\\) uniformément sur [0,1]</li>
              <li>calculer \\(Y = m(X) + \\varepsilon\\)</li>
              <li>générer \\(T\\) selon la loi spécifiée</li>
              <li>ne conserver que les triplets \\((X_i, Y_i, T_i)\\) tels que \\(Y_i \\geq T_i\\)</li>
            </ul>
            Elle renvoie à la fois l’échantillon complet et l’échantillon observé (après troncature).
          </p>
         
        <hr style='margin-top:40px;margin-bottom:20px;'>

        <h4><strong>📂 Code Source R Utilisé</strong></h4>
      ")),
        
        
        
        h5("1. Estimateur de Lynden-Bell \\( \\widehat{G}_{LB}(t) \\)"),
        pre(includeText("code/estimeLB_Gn.R")),
        
        h5("2. Estimation de \\( \\widehat{C}_n(y) \\)"),
        pre(includeText("code/C_n.R")),
        
        
        h5("4. Estimateur de Nadaraya-Watson corrigé"),
        pre(includeText("code/nw_estimator.R")),
        
        h5("5. Fonction de simulation des données tronquées"),
        pre(includeText("code/simulate_data.R")),
        
        
      )
    )
  ),
  
  nav_panel(
    "Estimation",
    card(
      card_header("Simulation Parameters", class = "bg-primary text-white"),
      layout_sidebar(
        sidebar = sidebar(
          sliderInput("n", "Sample size N:", min = 50, max = 1000, value = 300, step = 50),
          checkboxInput("auto_h", "Optimiser automatiquement h par validation croisée", value = FALSE),
          
          conditionalPanel(
            condition = "!input.auto_h",
            sliderInput("h", "Bandwidth h:", min = 0.01, max = 0.5, value = 0.1, step = 0.01)
          ),
          
          textInput("true_func", "True function m(x):", value = "sin(2 * pi * x)"),
          helpText("Example functions: x^2, sin(x), exp(x), log(x+1)"),
          
          # Méthode de troncature
          selectInput("trunc_method", "Méthode de troncature :",
                      choices = c("Mode aléatoire" = "random", 
                                  "Taux de troncature fixé" = "fixed"),
                      selected = "random"),
          
          # Panneau conditionnel pour taux fixe
          conditionalPanel(
            condition = "input.trunc_method == 'fixed'",
            sliderInput("trunc_rate", "Taux de troncature souhaité :", 
                        min = 0.1, max = 0.9, value = 0.5, step = 0.05)
          ),
          
          # Groupe des distributions uniquement visible en mode aléatoire
          conditionalPanel(
            condition = "input.trunc_method == 'random'",
            tagList(
              selectInput("t_dist", "Distribution de troncature :",
                          choices = c("Normale" = "normal",
                                      "Uniforme" = "uniform",
                                      "Exponentielle" = "exponential",
                                      "Gamma" = "gamma"),
                          selected = "normal"),
              
              conditionalPanel(
                condition = "input.t_dist == 'normal'",
                numericInput("t_mean", "Moyenne :", value = 0),
                numericInput("t_sd", "Écart-type :", value = 0.5, min = 0.1)
              ),
              
              conditionalPanel(
                condition = "input.t_dist == 'uniform'",
                numericInput("t_min", "Minimum :", value = -1),
                numericInput("t_max", "Maximum :", value = 1)
              ),
              
              conditionalPanel(
                condition = "input.t_dist == 'exponential'",
                numericInput("t_rate", "Taux :", value = 1, min = 0.1)
              ),
              
              conditionalPanel(
                condition = "input.t_dist == 'gamma'",
                numericInput("t_shape", "Forme :", value = 2, min = 0.1),
                numericInput("t_rate_gamma", "Taux :", value = 1, min = 0.1)
              )
            )
          ),
          
          numericInput("seed", "Random seed (optional):", value = NA),
          actionButton("simulate", "Run Simulation", class = "btn-success"),
          hr(),
          helpText("Adjust parameters and click the button to run a new simulation."),
          width = 300
        ),
        navset_card_tab(
          title = "Results",
          
          nav_panel(
            "Estimation Plot",
            plotlyOutput("main_plot", height = "500px"),
            verbatimTextOutput("selected_h_display")
            
          ),
          
          
          nav_panel(
            "Y Density (Before/After Truncation)",
            plotlyOutput("y_density_plot", height = "500px")
          ),
          
          nav_panel(
            "Data Summary",
            card(
              card_header("Truncation Summary"),
              tableOutput("trunc_summary"),
              card_header("Sample Data"),
              DTOutput("data_table")
            )
          )
        )
      )
    )
  ),nav_panel(
    "Performance",
    card(
      card_header("Bias, Variance and MSE Analysis", class = "bg-primary text-white"),
      layout_sidebar(
        sidebar = sidebar(
          sliderInput("n_sim", "Number of simulations:", min = 10, max = 200, value = 20),
          checkboxInput("auto_h_perf", "Optimiser h_n par validation croisée", value = TRUE),
          
          conditionalPanel(
            condition = "!input.auto_h_perf",
            sliderInput("h_perf", "Bande de lissage fixe (h):", min = 0.01, max = 0.5, value = 0.1)
          ),
          
          selectInput("n_values", "Sample sizes to compare:",
                      choices = c(50, 100, 200, 300, 500, 1000, 5000),
                      selected = c(100, 200, 500),
                      multiple = TRUE),
          actionButton("run_analysis", "Run Performance Analysis", class = "btn-success"),
          hr(),
          helpText("This analysis may take some time to complete depending on the number of simulations."),
          width = 300
        ),
        navset_card_tab(
          
          nav_panel(
            "Bias",
            plotlyOutput("bias_plot", height = "500px")
          ),
          nav_panel(
            "Variance",
            plotlyOutput("variance_plot", height = "500px")
          ),
          nav_panel(
            "MSE",
            plotlyOutput("mse_plot", height = "500px")
          ),
          
          
          nav_panel(
            "Summary",
            card(
              card_header("Performance Metrics Summary"),
              DTOutput("metrics_table"),
              card_header("Interpretation"),
              markdown("
              - **Bias**: Measures how far the average estimate is from the true value
              - **Variance**: Measures how much the estimates vary across simulations
              - **MSE**: Mean Squared Error (Bias² + Variance), overall accuracy measure
              
              A good estimator should have low bias, low variance, and thus low MSE.
              ")
            )
          )
        )
      )
    )
  ),
  nav_panel(
    "Bandwidth Study",
    card(
      card_header("Bandwidth Parameter (h) Study", class = "bg-primary text-white"),
      layout_sidebar(
        sidebar = sidebar(
          sliderInput("h_study_n", "Sample size:", min = 50, max = 1000, value = 300),
          sliderInput("h_range_study", "Bandwidth range:", 
                      min = 0.01, max = 1, value = c(0.03, 0.3), step = 0.01),
          sliderInput("h_steps", "Number of bandwidth values:", 
                      min = 3, max = 20, value = 10),
          actionButton("run_h_study", "Run Bandwidth Study", class = "btn-success"),
          hr(),
          helpText("This study shows how the bandwidth parameter affects the estimator."),
          width = 300
        ),
        navset_card_tab(
          title = "Results",
          nav_panel(
            "Multiple Plots",
            plotOutput("h_multiple_plots", height = "800px"),
            downloadButton("download_h_multiple", "Download PNG")
          ),
          nav_panel(
            "Interactive Plot",
            plotlyOutput("h_study_plot", height = "500px")
          ),
          nav_panel("MSE vs h", plotlyOutput("mse_vs_h_plot", height = "500px"),verbatimTextOutput("h_best_summary")),
          nav_panel("Bandwidth Summary Table", DT::DTOutput("bandwidth_table")),
          
          nav_panel(
            "Summary",
            card(
              card_header("Interpretation"),
              markdown("
              This panel studies how the bandwidth parameter h affects the estimator's performance.
              
              - **Small h**: Leads to more variance (overfitting)
              - **Large h**: Leads to more bias (underfitting)
              
              The optimal bandwidth balances bias and variance to minimize MSE.
              "),
              card_header("Optimal h"),
              
            )
          )
        )
      )
    )
  ),
  nav_panel(
    "MISE Simulation",
    card(
      card_header("MISE (Mean Integrated Squared Error) Simulation", class = "bg-primary text-white"),
      layout_sidebar(
        sidebar = sidebar(
          sliderInput("mise_n_sim", "Number of simulations:", min = 10, max = 200, value = 50),
          sliderInput("mise_n_values", "Sample sizes to test:", 
                      min = 50, max = 1000, value = c(100, 500), step = 50),
          actionButton("run_mise_simulation", "Run MISE Simulation", class = "btn-success"),
          hr(),
          helpText("This simulation may take some time. Lower number of simulations for faster results."),
          width = 300
        ),
        navset_card_tab(
          title = "Results",
          nav_panel(
            "Line Plots",
            plotlyOutput("mise_line_plot", height = "500px"),
            verbatimTextOutput("mise_best_summary")
          ),
          
          
          nav_panel("MISE Table", DT::DTOutput("mise_table")),
          
          
          
          nav_panel(
            "Summary",
            card(
              card_header("Interpretation"),
              markdown("
              This panel studies the Mean Integrated Squared Error (MISE) of the estimator.
              
              MISE measures the overall accuracy of the estimator by integrating the squared error over the domain.
              
              Key relationships to observe:
              - How MISE changes with bandwidth (h)
              - How MISE changes with sample size (n)
              
              Typically, we expect:
              - MISE decreases as n increases
              - MISE has a U-shaped relationship with h
              ")
            )
          )
        )
      )
    )
  ),
  
  nav_panel(
    "Confidence Intervals",
    card(
      card_header("Confidence Interval Simulation", class = "bg-primary text-white"),
      layout_sidebar(
        sidebar = sidebar(
          sliderInput("ci_n", "Sample size n:", min = 50, max = 2000, value = 500),
          
          checkboxInput("auto_h_ci", "Optimiser hₙ par validation croisée", value = TRUE),
          
          conditionalPanel(
            condition = "!input.auto_h_ci",
            sliderInput("ci_h", "Bande de lissage fixe (h):", min = 0.01, max = 0.5, value = 0.1)
          ),
          sliderInput("ci_b", "Number of simulations (B):", min = 10, max = 500, value = 100),
          sliderInput("ci_level", "Confidence level:", min = 0.8, max = 0.99, value = 0.95),
          actionButton("run_ci_simulation", "Run Simulation", class = "btn-success"),
          hr(),
          helpText("This simulation may take some time to complete. Lower number of simulations for faster results."),
          width = 300
        ),
        navset_card_tab(
          title = "Results",
          nav_panel(
            "Plot",
            plotlyOutput("ci_plot", height = "500px")
          ),
          nav_panel("Undercoverage Regions", plotlyOutput("ci_error_regions", height = "500px")),
          
          nav_panel(
            "Summary",
            card(
              card_header("Interpretation"),
              markdown("
            This panel simulates confidence intervals for the nonparametric estimator with truncation.
            
            - The dashed line shows the true function m(x)
            - The solid line shows the average of the estimates across simulations
            - The shaded region shows the empirical confidence intervals
            
            Good performance is indicated when:
            - The average estimate (solid line) closely follows the true function
            - The confidence band contains the true function across most of the domain
            ")
            )
          )
        )
      )
    )
  )
  
  
  
  
  
  
  ,
  
  nav_panel(
    "Uniform Convergence",
    card(
      card_header("Uniform Almost Sure Convergence",  class = "bg-primary text-white"),
      layout_sidebar(
        sidebar = sidebar(
          selectInput("unif_n_values", "Sample sizes to test:", 
                      choices = c(50, 100, 200, 500, 1000, 2000),
                      selected = c(100, 200, 500),
                      multiple = TRUE),
          checkboxInput("auto_h_unif", "Optimiser hₙ par validation croisée", value = TRUE),
          
          conditionalPanel(
            condition = "!input.auto_h_unif",
            sliderInput("unif_h", "Bande de lissage fixe (h):", min = 0.01, max = 0.5, value = 0.1)
          ),
          
          sliderInput("unif_n_sim", "Number of simulations:", min = 10, max = 500, value = 100),
          actionButton("run_unif_convergence", "Run Analysis", class = "btn-success"),
          hr(),
          helpText("This analysis may take some time for larger sample sizes and more simulations."),
          width = 300
        ),
        navset_card_tab(
          title = "Results",
          nav_panel(
            "Convergence Plot",
            plotlyOutput("unif_convergence_plot", height = "500px")
          ),
          nav_panel("Log-Log Fit", plotlyOutput("loglog_unif_plot", height = "500px")),
          
          
          
          nav_panel(
            "Error Matrix",
            DTOutput("unif_error_table")
          ),
          nav_panel(
            "Interpretation",
            card(
              card_header("Uniform Convergence Explanation"),
              markdown("
            This panel demonstrates uniform almost sure convergence of the estimator:
            
            - Plots the maximum absolute error (supremum norm) between the estimated and true functions
            - Shows how this error decreases as sample size increases
            - Uses multiple simulations to estimate the average behavior
            
            Key concepts:
            - **Uniform convergence**: The maximum error across all x values
            - **Almost sure convergence**: The error converges to 0 with probability 1
            - The log-scale x-axis helps visualize the rate of convergence
            
            In well-behaved cases, we expect to see the error decreasing as sample size increases.
            ")
            )
          )
        )
      )
    )
  ),
  nav_panel(
    "Estimator Comparison",
    card(
      card_header("Comparaison entre estimateurs non paramétriques", class = "bg-primary text-white"),
      navset_card_tab(
        title = "Sections",
        
        nav_panel(
          "Théorie",
          card(
            card_header("Méthodologie : Estimateurs Non Paramétriques", class = "bg-primary text-white"),
            card_body(
              withMathJax(HTML("
<h4><strong>📚 Estimateur Spline</strong></h4>
<p>Les splines ajustent une fonction \\( m(x) \\) en reliant des morceaux de polynômes avec continuité des dérivées. Le modèle est :</p>
\\[
\\widehat{m}(x) = \\arg\\min_{f} \\left\\{ \\sum_{i=1}^n (Y_i - f(X_i))^2 + \\lambda \\int (f''(x))^2 dx \\right\\}
\\]
<p>où \\( \\lambda \\) est le paramètre de lissage (lié à <code>spar</code> dans R).</p>

<h4><strong>📚 Estimateur LOESS (Local Regression)</strong></h4>
<p>LOESS effectue une régression locale pondérée autour de chaque point :</p>
\\[
\\widehat{m}(x) = \\arg\\min_{\\beta_0, \\beta_1} \\sum_{i=1}^n K_h(x - X_i) (Y_i - \\beta_0 - \\beta_1 (X_i - x))^2
\\]
<p>où le noyau \\( K_h \\) donne plus de poids aux points proches. Le paramètre <code>span</code> contrôle la fenêtre de lissage.</p>

<h4><strong>📚 Estimateur Polynomial Local</strong></h4>
<p>L'estimateur polynomial local d'ordre \\( p \\) ajuste un polynôme autour de chaque point :</p>
\\[
\\widehat{m}(x) = \\mathbf{e}_1^T (X^T W X)^{-1} X^T W \\mathbf{Y}
\\]
<p>où :</p>
<ul>
<li>\\( X \\) est la matrice des puissances centrées de \\( X_i \\)</li>
<li>\\( W \\) est une matrice diagonale de poids définis par un noyau \\( K_h(x - X_i) \\)</li>
<li>\\( \\mathbf{e}_1 \\) sélectionne l’ordonnée à l’origine du polynôme local</li>
</ul>

<hr>
<h4><strong>📂 Code Source R Utilisé</strong></h4>
<pre style='background:#f8f9fa;border:1px solid #dee2e6;padding:10px;overflow-x:auto'>
<code>", paste(readLines("code/np_estimator_comparison.R"), collapse = "\n"), "</code>
</pre>
      "))
            )
          )
        ),
        nav_panel(
          "Plots",
          
          card(
            card_header("Comparison Nonparametric Estimators", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                sliderInput("comp_n", "Sample size n:", min = 50, max = 2000, value = 500),
                checkboxInput("auto_h_comp", "Optimiser hₙ pour NW (CV)", value = TRUE),
                conditionalPanel(
                  condition = "!input.auto_h_comp",
                  sliderInput("comp_h", "Bandwidth h (for NW):", min = 0.01, max = 0.5, value = 0.1)
                ),
                sliderInput("comp_loess_span", "LOESS span:", min = 0.1, max = 1, value = 0.3),
                sliderInput("comp_poly_deg", "Local polynomial degree:", min = 1, max = 3, value = 2),
                sliderInput("comp_spline_spar", "Spline smoothing parameter:", min = 0.1, max = 1, value = 0.7),
                
                actionButton("run_comparison", "Run Comparison", class = "btn-success"),
                hr(),
                helpText("This comparison includes MSE and relative efficiency metrics for multiple estimators."),
                width = 300
              ),
              plotOutput("comparison_plot", height = "800px"),
              # Dans l'UI
              downloadButton("download_comparison_plot", "Télécharger la figure")
              
            )
          )
          
        ),
        nav_panel(
          "MSE Comparison",
          card(
            card_header("MSE barplot", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                sliderInput("comp_n", "Sample size n:", min = 50, max = 2000, value = 500),
                checkboxInput("auto_h_comp", "Optimiser hₙ pour NW (CV)", value = TRUE),
                conditionalPanel(
                  condition = "!input.auto_h_comp",
                  sliderInput("comp_h", "Bandwidth h (for NW):", min = 0.01, max = 0.5, value = 0.1)
                ),
                sliderInput("comp_loess_span", "LOESS span:", min = 0.1, max = 1, value = 0.3),
                sliderInput("comp_poly_deg", "Local polynomial degree:", min = 1, max = 3, value = 2),
                sliderInput("comp_spline_spar", "Spline smoothing parameter:", min = 0.1, max = 1, value = 0.7),
                
                actionButton("run_comparison", "Run Comparison", class = "btn-success"),
                hr(),
                helpText("This comparison includes MSE and relative efficiency metrics for multiple estimators."),
                width = 300
              ),
              plotlyOutput("mse_comparison_plot", height = "500px")
            )
            
          )
        ),
        
        
        nav_panel(
          "Interpretation", 
          card(
            card_header("Interpretation", class = "bg-primary text-white"),
            markdown(
              "
### 📌 Estimators Compared:
- **Nadaraya-Watson (corrected)**: Kernel-based estimator that corrects for truncation bias
- **LOESS**: Local polynomial regression with weighted smoothing
- **Local Polynomial**: Polynomial fit of user-defined degree
- **Splines**: Spline smoothing using regularization

### 📈 Metrics:
- **MSE**: Mean Squared Error over the grid




### 🔍 Recommendations:

- Be mindful of **overfitting** with low bandwidth or high degree
- Consider **spline or local polynomial** for smoother or more structured underlying functions
            ")
          ))
        
        
      )
      
      
    )
  ),
  
  
  
  # Ajouter dans la navbar (après le dernier nav_panel)
  
  
  nav_panel(
    "Simple NN + NW Mixture",
    card(
      card_header("Combinaison des Estimateurs Nadaraya-Watson & Réseau de Neurones", class = "bg-primary text-white"),
      
      navset_card_tab(
        title = "Sections",
        
        # Onglet 1 : Explication
        nav_panel(
          "Théorie et Code R",
          card_header("Méthodologie & Code R", class = "bg-primary text-white"),
          card_body(
            withMathJax(HTML("
            
<h4><strong>🧮 Formulation mathématique</strong></h4>

<p>Un réseau de neurones à une couche cachée et une sortie linéaire (cas de la régression) se définit comme suit :</p>

<ul>
  <li>Entrée : \\( \\mathbf{X} \\in \\mathbb{R}^d \\)</li>
  <li>Couche cachée : \\( H = \\texttt{size} \\) neurones, activation \\( \\sigma(u) = \\frac{1}{1 + e^{-u}} \\)</li>
  <li>Sortie : estimation \\( \\widehat{Y} \\) corrigée pour troncature</li>
</ul>

<p>La fonction estimée est :</p>

\\[
\\widehat{Y} = \\sum_{j=1}^H w_j^{(2)} \\cdot \\sigma\\left( \\sum_{i=1}^d w_{ij}^{(1)} X_i + b_j^{(1)} \\right) + b^{(2)}
\\]

<ul>
  <li>\\( w_{ij}^{(1)} \\) : poids entrée vers couche cachée</li>
  <li>\\( b_j^{(1)} \\) : biais des neurones cachés</li>
  <li>\\( w_j^{(2)} \\) : poids couche cachée vers sortie</li>
  <li>\\( b^{(2)} \\) : biais de sortie</li>
</ul>

<p>La fonction de coût avec régularisation (weight decay) devient :</p>

\\[
L = \\sum_{i=1}^n w_i \\cdot (Y_i - \\widehat{Y}_i)^2 + \\texttt{decay} \\cdot \\left( \\|W^{(1)}\\|_2^2 + \\|W^{(2)}\\|_2^2 \\right)
\\]

<p>où \\( w_i = \\frac{1}{\\widehat{G}_{LB}(Y_i)} \\) corrige la troncature à gauche.</p>

<hr>



<h4><strong>🧠 Schématisation du Réseau de Neurones</strong></h4>

<div style='display:flex;justify-content:space-between;align-items:flex-start;gap:40px;'>

  <div style='flex:1.75; text-align:center; '>
    <img src='réseau.png' style='width:80%; max-width:2500px; border:1px solid #ccc; padding:10px; background:#fff'/>
  </div>

  <div style='flex:2; align-self: center; margin-left: 130px;'>
  <p><strong>Partie haute – Architecture</strong></p>
  <ul>
    <li>Entrée : \\( X_1, ..., X_d \\) vers chaque neurone caché</li>
    <li>Couche cachée : \\( \\texttt{size} \\) neurones avec activation sigmoïde</li>
    <li>Sortie : combinaison linéaire pour estimer \\( \\widehat{Y} \\)</li>
    <li>Troncature : seules les observations avec \\( Y_i \\geq T_i \\) sont conservées</li>
  </ul>

  <p><strong>Partie basse – Pipeline d'entraînement</strong></p>
  <ol>
    <li>Simulation des données tronquées \\( (X_i, Y_i, T_i) \\)</li>
    <li>Estimation \\( \\widehat{G}_{LB}(Y_i) \\)</li>
    <li>Poids \\( w_i = 1 / \\widehat{G}_{LB}(Y_i) \\)</li>
    <li>Apprentissage via <code>nnet()</code> avec <code>size</code>, <code>decay</code>, <code>maxit</code></li>
    <li>Prédiction sur une grille \\( x \\in [0,1] \\) pour estimer \\( m(x) \\)</li>
  </ol>
</div>

</div>
<p>Cette illustration montre comment intégrer la correction de troncature dans un réseau de neurones à une couche cachée, entraîné en régression via le package <code>nnet</code>.</p>
<p>La souplesse des réseaux est ainsi combinée à la robustesse de la pondération via l’estimateur de Lynden-Bell.</p>


<hr>


<h4><strong>🔧 Paramètres utilisés dans l'estimation</strong></h4>

<table border='1' style='border-collapse: collapse;'>
  <thead>
    <tr><th>Paramètre</th><th>Interprétation</th></tr>
  </thead>
  <tbody>
    <tr><td><code>size</code></td><td>Nombre de neurones dans la couche cachée</td></tr>
    <tr><td><code>decay</code></td><td>Régularisation L2 pour éviter le surapprentissage</td></tr>
    <tr><td><code>maxit</code></td><td>Nombre d'itérations pour l'entraînement</td></tr>
  
  </tbody>
</table>

<hr>

<h4><strong>🧠 Méthode de Combinaison NN + NW</strong></h4>

<p>
Nous combinons deux estimateurs de la fonction de régression \\( m(x) = \\mathbb{E}[Y \\mid X = x] \\) dans le cadre de données tronquées à gauche :
</p>

<ul>
  <li>\\( \\widehat{m}_{NW}(x) \\) : estimateur de Nadaraya-Watson corrigé pour troncature</li>
  <li>\\( \\widehat{m}_{NN}(x) \\) : estimateur basé sur un réseau de neurones à une couche cachée</li>
</ul>

<p>
L'estimateur final est une <strong>combinaison convexe</strong> des deux :
</p>

<p>
\\[
\\widehat{m}_{\\text{mix}}(x) = \\lambda \\widehat{m}_{NW}(x) + (1 - \\lambda) \\widehat{m}_{NN}(x), \\quad \\lambda \\in [0, 1]
\\]
</p>

<p>
Le paramètre \\( \\lambda \\) est ajusté par l'utilisateur pour pondérer l'importance des deux méthodes.
</p>

<hr>



<h4><strong>📂 Code Source R Utilisé</strong></h4>

<h5>1. Simulation de données tronquées</h5>
<pre style='background:#f8f9fa;border:1px solid #dee2e6;padding:10px;'>
<code>
", paste(readLines("code/simulate_data.R"), collapse = "\n"), "
</code>
</pre>

<h5>2. Estimateur Nadaraya-Watson corrigé</h5>
<pre style='background:#f8f9fa;border:1px solid #dee2e6;padding:10px;'>
<code>
", paste(readLines("code/nw_estimator.R"), collapse = "\n"), "
</code>
</pre>

<h5>3. Entraînement du Réseau de Neurones</h5>
<pre style='background:#f8f9fa;border:1px solid #dee2e6;padding:10px;'>
<code>
", paste(readLines("code/train_nn_estimator.R"), collapse = "\n"), "
</code>
</pre>

          "))
          )
        ),
        nav_panel(
          "Comparaison séparée",
          card(
            card_header("Comparaison des Estimateurs Séparée", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                sliderInput("mix_lambda", "Poids lambda :", min = 0, max = 1, value = 0.5, step = 0.05),
                checkboxInput("auto_h_mix", "Optimiser hₓ par validation croisée", value = TRUE),
                
                conditionalPanel(
                  condition = "!input.auto_h_mix",
                  sliderInput("mix_h", "Bande de lissage h:", min = 0.01, max = 0.5, value = 0.1, step = 0.01)
                ),
                numericInput("nn_sample_size", "Taille de l'échantillon n:", value = 300, min = 50, max = 2000, step = 50),
                sliderInput("nn_size", "Nombre de neurones cachés:", min = 1, max = 50, value = 10),
                sliderInput("nn_decay", "Décroissance (weight decay):", min = 0, max = 0.1, value = 0.01, step = 0.001),
                sliderInput("nn_maxit", "Nombre maximal d'itérations:", min = 100, max = 5000, value = 1000),
                sliderInput("nn_seed", "Seed (réplicabilité):", min = 1, max = 9999, value = 123),
                actionButton("run_mix_plot", "Lancer la combinaison", class = "btn-success"),
                width=300
              ),
              plotlyOutput("nn_mix_comparison_plot", height = "600px")
            )
          )
        ),
        nav_panel(
          "Étude Lambda",
          card(
            card_header("Étude du paramètre de combinaison \u03bb", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                numericInput("lambda_n", "Taille échantillon (n):", value = 300, min = 50, max = 2000, step = 50),
                checkboxInput("auto_h_lambda", "Optimiser hₙ par validation croisée", value = TRUE),
                
                conditionalPanel(
                  condition = "!input.auto_h_lambda",
                  sliderInput("lambda_h", "Bande de lissage h:", min = 0.01, max = 0.5, value = 0.1)
                ),
                
                sliderInput("lambda_nn_size", "Neurones cachés:", min = 1, max = 50, value = 10),
                sliderInput("lambda_nn_decay", "Weight decay:", min = 0, max = 0.1, value = 0.01, step = 0.001),
                sliderInput("lambda_nn_maxit", "Max iterations:", min = 100, max = 5000, value = 1000),
                sliderInput("lambda_nsim", "Nombre de simulations:", min = 5, max = 100, value = 20),
                actionButton("run_lambda_study", "Lancer l'étude", class = "btn-success"),
                hr()
                
              ),
              plotlyOutput("lambda_mse_plot", height = "500px"),
              verbatimTextOutput("best_lambda_display")
            )
          )
        ),
        nav_panel(
          "Paramètres NN",
          card(
            card_header("Étude des paramètres du réseau de neurones", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                selectInput("nn_param_to_study", "Paramètre à étudier:", choices = c("size", "decay", "maxit"), selected = "size"),
                sliderInput("nn_study_lambda", "Lambda (\u03bb):", min = 0, max = 1, value = 0.5, step = 0.05),
                checkboxInput("auto_h_nn_param", "Optimiser hₙ par validation croisée", value = TRUE),
                
                conditionalPanel(
                  condition = "!input.auto_h_nn_param",
                  sliderInput("nn_study_h", "Bande de lissage h:", min = 0.01, max = 0.5, value = 0.1, step = 0.01)
                ),
                
                numericInput("nn_study_n", "Taille de l'échantillon n:", value = 300, min = 50, max = 2000, step = 50),
                sliderInput("nn_study_nsim", "Nombre de simulations:", min = 5, max = 100, value = 20),
                sliderInput("nn_fixed_size", "Taille couche cachée (fixe):", min = 1, max = 50, value = 10),
                sliderInput("nn_fixed_decay", "Weight decay (fixe):", min = 0, max = 0.1, value = 0.01, step = 0.001),
                sliderInput("nn_fixed_maxit", "Max iterations (fixe):", min = 100, max = 5000, value = 1000),
                actionButton("run_nn_param_study", "Lancer l'étude", class = "btn-success"),
                width = 300,
                hr()
                
              ),
              plotlyOutput("nn_param_study_plot", height = "600px"),
              verbatimTextOutput("best_nn_param_display")
            )
          )
        ),
        nav_panel(
          "Paramètres optimaux",
          card(
            card_header("Résumé des paramètres optimaux pour la combinaison NN + NW", class = "bg-primary text-white"),
            card_body(
              DT::dataTableOutput("summary_optimal_params"),
              br(),
              verbatimTextOutput("summary_optimal_params_text")
            )
          )
        ),
        nav_panel(
          "NN seul - tuning",
          card(
            card_header("Étude des Paramètres du Réseau de Neurones Seul", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                selectInput("nn_solo_param", "Paramètre à étudier:", choices = c("size", "decay", "maxit"), selected = "size"),
                numericInput("nn_solo_n", "Taille de l'échantillon:", value = 300, min = 50, max = 2000, step = 50),
                sliderInput("nn_solo_nsim", "Nombre de simulations:", min = 5, max = 100, value = 20),
                sliderInput("nn_solo_fixed_size", "Taille couche cachée (fixe):", min = 1, max = 50, value = 10),
                sliderInput("nn_solo_fixed_decay", "Weight decay (fixe):", min = 0, max = 0.1, value = 0.01, step = 0.001),
                sliderInput("nn_solo_fixed_maxit", "Max iterations (fixe):", min = 100, max = 5000, value = 1000),
                actionButton("run_nn_solo_study", "Lancer l'étude", class = "btn-success"),
                width = 300,
                hr()
              ),
              plotlyOutput("nn_solo_param_plot", height = "600px"),
              verbatimTextOutput("best_nn_solo_param_display")
            )
          )
        ),
        
        # --- Résumé tabulaire ---
        nav_panel(
          "Récapitulatif NN seul",
          card(
            card_header("Paramètres optimaux du réseau seul", class = "bg-primary text-white"),
            card_body(
              DT::dataTableOutput("summary_nn_solo_params")
            )
          )
        ),
        nav_panel(
          "Comparaison finale",
          card(
            card_header("Comparaison des Estimateurs NW, NN et Mixture", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                sliderInput("final_h_nw", "h pour NW:", min = 0.01, max = 0.5, value = 0.1, step = 0.01),
                sliderInput("final_size_nn", "Neurones cachés NN seul:", min = 1, max = 50, value = 10),
                sliderInput("final_decay_nn", "Decay NN seul:", min = 0, max = 0.1, value = 0.01, step = 0.001),
                sliderInput("final_maxit_nn", "Maxit NN seul:", min = 100, max = 5000, value = 1000),
                
                hr(),
                sliderInput("final_lambda_mix", "\u03bb (mixture):", min = 0, max = 1, value = 0.5, step = 0.05),
                sliderInput("final_h_mix", "h (mixture):", min = 0.01, max = 0.5, value = 0.1, step = 0.01),
                sliderInput("final_size_nn_mix", "Neurones cachés (mixture):", min = 1, max = 50, value = 10),
                sliderInput("final_decay_nn_mix", "Decay (mixture):", min = 0, max = 0.1, value = 0.01, step = 0.001),
                sliderInput("final_maxit_nn_mix", "Maxit (mixture):", min = 100, max = 5000, value = 1000),
                
                numericInput("final_n", "Taille de l\u2019\u00e9chantillon:", value = 300, min = 50, max = 2000, step = 50),
                sliderInput("final_nsim", "Nombre de simulations:", min = 5, max = 100, value = 20),
                
                actionButton("run_final_comparison", "Lancer la comparaison", class = "btn-success"),
                width = 300
              ),
              layout_column_wrap(
                width = 1,
                plotlyOutput("final_comparison_plot", height = "500px"),
                DT::dataTableOutput("final_comparison_table")
              )
            )
          )
        )
        
      )
    )
    
  ),
  nav_panel(
    "Méthode de Stacking / Apprentissage Méta",
    card(
      card_header("Combinaison des Estimateurs", class = "bg-primary text-white"),
      
      navset_card_tab(
        title = "Sections",
        
        nav_panel(
          "Théorie et code",
          card(
            card_header("Méthodologie & Code R", class = "bg-primary text-white"),
            card_body(
              withMathJax(HTML("
<h4><strong>🔀 Méthode de Stacking / Apprentissage Méta</strong></h4>

<div style='font-size: 1.05em;'>

<h5>🧠 Idée générale</h5>
<p>
L’approche de <strong>stacking</strong> vise à <em>combiner intelligemment plusieurs estimateurs</em> \\( \\widehat{m}^{(1)}(x), \\widehat{m}^{(2)}(x), ..., \\widehat{m}^{(K)}(x) \\) 
via un modèle superviseur (appelé <em>meta-learner</em>) qui apprend la meilleure combinaison en minimisant l’erreur sur un jeu de validation.
</p>

<h5>⚙️ Dans notre cas</h5>
<p>On souhaite apprendre une fonction \\( f_{\\text{meta}} \\) telle que :</p>
\\[
\\widehat{m}_{\\text{meta}}(x) = f_{\\text{meta}}\\left( \\widehat{m}_{NW}(x), \\widehat{m}_{NN}(x) \\right)
\\]
<p>où :</p>
<ul>
  <li>\\( \\widehat{m}_{NW}(x) \\) est l’estimation par <strong>Nadaraya-Watson corrigé</strong></li>
  <li>\\( \\widehat{m}_{NN}(x) \\) est l’estimation par <strong>réseau de neurones corrigé</strong></li>
  <li>\\( f_{\\text{meta}} \\) est un modèle à apprendre (régression linéaire, réseau simple, etc.)</li>
</ul>

<h5>📐 Étapes de la méthode</h5>
<ol>
  <li><strong>Division des données</strong><br>
      On divise les données observées \\( (X_i, Y_i) \\) avec \\( Y_i \\geq T_i \\) en :
      <ul>
        <li>Un ensemble d’apprentissage : \\( \\mathcal{D}_{\\text{train}} \\)</li>
        <li>Un ensemble de validation : \\( \\mathcal{D}_{\\text{val}} \\)</li>
      </ul>
  </li>

  <li><strong>Apprentissage des estimateurs de base sur \\( \\mathcal{D}_{\\text{train}} \\)</strong><br>
      On ajuste :
      \\[
      \\widehat{m}^{\\text{train}}_{NW}(x), \\quad \\widehat{m}^{\\text{train}}_{NN}(x)
      \\]
  </li>

  <li><strong>Construction des prédicteurs sur \\( \\mathcal{D}_{\\text{val}} \\)</strong><br>
      Pour chaque point \\( x_j \\in \\mathcal{D}_{\\text{val}} \\), on calcule :
      \\[
      \\widehat{m}_{NW,j}, \\quad \\widehat{m}_{NN,j}
      \\]
  </li>

  <li><strong>Apprentissage du modèle méta</strong><br>
      À partir des paires \\( ((\\widehat{m}_{NW,j}, \\widehat{m}_{NN,j}), Y_j) \\), on ajuste :
      \\[
      f_{\\text{meta}}(a, b) \\approx Y_j
      \\]
      <p>Par exemple, si \\( f_{\\text{meta}} \\) est une régression linéaire :</p>
      \\[
      \\widehat{m}_{\\text{meta}}(x) = \\beta_0 + \\beta_1 \\widehat{m}_{NW}(x) + \\beta_2 \\widehat{m}_{NN}(x)
      \\]
  </li>

  <li><strong>Prédiction finale</strong><br>
      Sur une nouvelle grille de \\( x \\), on prédit :
      \\[
      \\widehat{m}_{\\text{meta}}(x) = f_{\\text{meta}}\\left( \\widehat{m}^{\\text{train}}_{NW}(x), \\widehat{m}^{\\text{train}}_{NN}(x) \\right)
      \\]
  </li>
</ol>

<h5>✅ Avantages du stacking</h5>
<ul>
  <li>Combine automatiquement les forces des estimateurs de base (ex: biais faible de NW, flexibilité du NN)</li>
  <li>Souvent plus performant que chaque méthode prise individuellement</li>
  <li>Permet d’utiliser des modèles non-linéaires pour \\( f_{\\text{meta}} \\) (réseaux de neurones, arbres, etc.)</li>
</ul>

</div>

<h5>📚 <strong>Rappels mathématiques des modèles de stacking :</strong></h5>

<p><strong>Régression Linéaire :</strong><br>
L'estimateur de Régression Linéaire est donné par :</p>
\\[
\\widehat{Y} = \\beta_0 + \\sum_{j=1}^{p} \\beta_j X_j
\\]

<p><strong>Support Vector Regression (SVR) :</strong><br>
L'estimateur SVR est donné par :</p>
\\[
\\min_{w, b, \\xi, \\xi^*} \\frac{1}{2} \\|w\\|^2 + C \\sum_{i=1}^{n} (\\xi_i + \\xi_i^*)
\\]
<p>sous contraintes :</p>
\\[
\\begin{aligned}
& Y_i - w^T \\phi(X_i) - b \\le \\epsilon + \\xi_i \\\\
& w^T \\phi(X_i) + b - Y_i \\le \\epsilon + \\xi_i^* \\\\
& \\xi_i, \\xi_i^* \\ge 0
\\end{aligned}
\\]

<p><strong>XGBoost :</strong><br>
L'estimateur XGBoost est donné par :</p>
\\[
\\text{Obj}(t) = \\sum_{i=1}^n l(y_i, \\hat{y}_i^{(t-1)} + f_t(x_i)) + \\Omega(f_t)
\\]
<p>avec \\( \\Omega(f) = \\gamma T + \\frac{1}{2} \\lambda \\|w\\|^2 \\), où \\( f_t \\) est un arbre, \\( l \\) la fonction de perte.</p>

<p><strong>MARS (Multivariate Adaptive Regression Splines) :</strong><br>
L'estimateur MARS est donné par :</p>
\\[
\\widehat{Y} = \\beta_0 + \\sum_{m=1}^M \\beta_m B_m(X)
\\]
<p>où \\( B_m(X) \\) sont des fonctions de base (splines adaptatifs) et \\( \\beta_m \\) les coefficients.</p>

<p><strong>Random Forest :</strong><br>
L'estimateur Random Forest est donné par :</p>
\\[
\\widehat{Y} = \\frac{1}{T} \\sum_{t=1}^{T} f_t(X)
\\]
<p>où \\( f_t \\) est un arbre de décision entraîné sur un sous-échantillon bootstrap de l'ensemble des données.</p>
<h4><strong>📂 Code Source R Utilisé</strong></h4>

")),
              h5("1. Simulation et Stacking des estimateurs \\( \\widehat{m}_{NW}, \\widehat{m}_{NN} \\) et apprentissage méta"),
              pre(includeText("Code/stacking.R"))
              
            )
          )
        )
        ,
        
        nav_panel(
          "Visualisation Stacking Séparée",
          card(
            card_header("Stacking: Comparaison des Estimateurs", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                numericInput("stack_n", "Taille de l'échantillon:", value = 300, min = 50, max = 2000, step = 50),
                
                checkboxInput("auto_h_stack", "Optimiser hₙ par validation croisée", value = TRUE),
                conditionalPanel(
                  
                  condition = "!input.auto_h_stack",
                  sliderInput("stack_h", "Bande de lissage h (NW):", min = 0.01, max = 0.5, value = 0.1)
                ),
                
                # Paramètres NN
                sliderInput("stack_nn_size", "Neurones cachés NN:", min = 1, max = 50, value = 10),
                sliderInput("stack_nn_decay", "Weight decay:", min = 0, max = 0.1, value = 0.01, step = 0.001),
                sliderInput("stack_nn_maxit", "Max iterations NN:", min = 100, max = 5000, value = 1000),
                
                # Sélection des méthodes pour le stacking
                selectInput("stacking_method", "Méthode stacking meta-learner:",
                            choices = c("Régression linéaire" = "lm", 
                                        "MARS" = "mars", 
                                        "SVR" = "svr", 
                                        "XGBoost" = "xgb",
                                        "Random Forest" = "rf"),
                            selected = "lm"),
                
                checkboxGroupInput("stack_methods", "Méthodes à inclure dans le stacking:",
                                   choices = c("Nadaraya-Watson" = "nw",
                                               "Réseau de neurones" = "nn"),
                                   selected = c("nw", "nn")),
                
                actionButton("run_stack_plot", "Lancer la simulation", class = "btn-success"),
                width = 300
              ),
              plotlyOutput("stack_plot_facet")
            )
          )
        )
        ,
        
        nav_panel(
          "Optimisation Paramètres du Stacking",
          card(
            card_header("Optimisation des Paramètres du Stacking", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                selectInput("stack_param", "Paramètre à étudier:", choices = c("size", "decay", "maxit"), selected = "size"),
                sliderInput("stack_nsim", "Nombre de simulations:", min = 5, max = 100, value = 20),
                
                numericInput("stack_n", "Taille de l'échantillon:", value = 300, min = 50, max = 2000, step = 50),
                checkboxInput("auto_h_stack", "Optimiser hₙ par validation croisée", value = TRUE),
                conditionalPanel(condition = "!input.auto_h_stack",
                                 sliderInput("stack_h", "Bande de lissage h (NW):", min = 0.01, max = 0.5, value = 0.1)),
                
                # Paramètres fixes
                sliderInput("stack_fixed_size", "Taille couche cachée NN:", min = 1, max = 50, value = 10),
                sliderInput("stack_fixed_decay", "Weight decay NN:", min = 0, max = 0.1, value = 0.01, step = 0.001),
                sliderInput("stack_fixed_maxit", "Max iterations NN:", min = 100, max = 5000, value = 1000),
                
                # Sélection stacking
                selectInput("stacking_method", "Méthode stacking meta-learner:",
                            choices = c("lm", "mars", "svr", "xgb", "rf"), selected = "lm"),
                checkboxGroupInput("stack_methods", "Méthodes à inclure dans le stacking:",
                                   choices = c("Nadaraya-Watson" = "nw", "Réseau de neurones" = "nn"),
                                   selected = c("nw", "nn")),
                
                actionButton("run_stack_tuning", "Lancer l'optimisation", class = "btn-success"),
                hr()
                
              ),
              plotlyOutput("stack_param_plot"),
              verbatimTextOutput("best_stack_param_display")
            )
          )
        )
        ,
        
        nav_panel(
          "Récapitulatif Stacking",
          card(
            card_header("Paramètres optimaux du modèle Stacking", class = "bg-primary text-white"),
            card_body(
              DT::dataTableOutput("summary_stacking_params"),
              br(),
              verbatimTextOutput("summary_stacking_text")
            )
          )
        ),
        
        nav_panel(
          "Tuning NN seul (Stacking)",
          card(
            card_header("Tuning du Réseau NN utilisé dans le Stacking", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                selectInput("stacking_nn_param", "Paramètre à étudier:", choices = c("size", "decay", "maxit")),
                numericInput("stacking_nn_n", "Taille de l’échantillon:", value = 300, min = 50, max = 2000, step = 50),
                sliderInput("stacking_nn_nsim", "Nombre de simulations:", min = 5, max = 100, value = 20),
                sliderInput("stacking_nn_fixed_size", "Taille (fixe):", min = 1, max = 50, value = 10),
                sliderInput("stacking_nn_fixed_decay", "Decay (fixe):", min = 0, max = 0.1, value = 0.01),
                sliderInput("stacking_nn_fixed_maxit", "Maxit (fixe):", min = 100, max = 5000, value = 1000),
                actionButton("run_stacking_nn_study", "Lancer", class = "btn-success"),
                hr()
              ),
              plotlyOutput("stacking_nn_param_plot", height = "600px"),
              verbatimTextOutput("best_stacking_nn_param_display")
            )
          )
        ),
        
        nav_panel(
          "Récapitulatif NN (Stacking)",
          card(
            card_header("Résumé NN Stacking", class = "bg-primary text-white"),
            card_body(
              DT::dataTableOutput("summary_stacking_nn_params")
            )
          )
        ),
        nav_panel(
          "Comparaison finale Stacking",
          card(
            card_header("Comparaison finale des méthodes", class = "bg-primary text-white"),
            layout_sidebar(
              sidebar = sidebar(
                numericInput("final_stack_n", "Taille de l'échantillon:", value = 300, min = 50, max = 2000, step = 50),
                sliderInput("final_stack_nsim", "Nombre de simulations:", min = 5, max = 100, value = 20),
                
                checkboxInput("final_auto_h", "Optimiser hₙ (NW) par validation croisée", value = TRUE),
                conditionalPanel("!input.final_auto_h", sliderInput("final_stack_h", "Bande lissage h:", 0.01, 0.5, 0.1)),
                
                h4("Paramètres NN seul"),
                sliderInput("final_nn_size_solo", "Taille NN", 1, 50, 10),
                sliderInput("final_nn_decay_solo", "Decay NN", 0, 0.1, 0.01),
                sliderInput("final_nn_maxit_solo", "Maxit NN", 100, 5000, 1000),
                
                h4("Paramètres NN dans Stacking"),
                sliderInput("final_nn_size_stack", "Taille NN", 1, 50, 10),
                sliderInput("final_nn_decay_stack", "Decay NN", 0, 0.1, 0.01),
                sliderInput("final_nn_maxit_stack", "Maxit NN", 100, 5000, 1000),
                
                h4("XGBoost et RF"),
                sliderInput("final_rf_ntree", "RF: ntree", 50, 1000, 500),
                sliderInput("final_xgb_nrounds", "XGBoost: nrounds", 50, 1000, 100),
                
                selectInput("stacking_method", "Méthode stacking meta-learner:",
                            choices = c("lm", "mars", "svr", "xgb", "rf"),
                            selected = "lm"),
                
                checkboxGroupInput("stack_methods", "Méthodes à inclure dans le stacking :",
                                   choices = c("Nadaraya-Watson" = "nw", "Réseau de neurones" = "nn"),
                                   selected = c("nw", "nn")),
                
                actionButton("run_final_stack_comparison", "Lancer", class = "btn-success"),
                width = 300
              ),
              layout_column_wrap(width = 1,
                                 plotlyOutput("final_stack_plot", height = "600px"),
                                 DT::dataTableOutput("final_stack_table"))
            )
          )
        )
        
        
        
        
      )
    )
  )
  
  
  
  
)


# Server
server <- function(input, output, session) {
  
  # Parse the true function from user input
  parse_true_func <- reactive({
    function(x) {
      tryCatch({
        eval(parse(text = input$true_func), envir = list(x = x))
      }, error = function(e) {
        # Default to sine function if parsing fails
        sin(2 * pi * x)
      })
    }
  })
  
  # Get truncation distribution parameters
  get_t_params <- reactive({
    if (input$trunc_method == "fixed") {
      return(c(input$trunc_rate))
    } else {
      switch(input$t_dist,
             "normal" = c(input$t_mean, input$t_sd),
             "uniform" = c(input$t_min, input$t_max),
             "exponential" = c(input$t_rate),
             "gamma" = c(input$t_shape, input$t_rate_gamma),
             c(0, 0.5)) # valeur par défaut
    }
  })
  
  
  sim_data <- eventReactive(input$simulate, {
    seed <- if (is.na(input$seed)) NULL else input$seed
    m_true <- parse_true_func()
    T_params <- get_t_params()
    simulate_data(input$n, m_true, input$t_dist, T_params, input$trunc_method, seed)
  })
  
  est_results <- reactive({
    withProgress(message = 'Calcul de l\'estimateur en cours...', value = 0, {
      
      data <- sim_data()
      m_true <- parse_true_func()
      x_grid <- seq(0, 1, length.out = 100)
      
      incProgress(0.3, detail = "Choix de la bande de lissage h")
      # Bande automatique si activée
      if (input$auto_h) {
        h_grid <- seq(0.03, 0.4, length.out = 30)
        h_val <- select_h_cv(data$X_obs, data$Y_obs, data$T_obs, h_grid)
      } else {
        h_val <- input$h
      }
      
      incProgress(0.6, detail = "Estimation de la fonction m(x)")
      m_est <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h_val)
      
      incProgress(1, detail = "Finalisation")
      list(
        x = x_grid,
        m_est = m_est,
        m_true = m_true(x_grid),
        data = data,
        h_val = h_val
      )
    })
  })
  
  
  
  optimal_params <- reactiveValues(
    lambda = NA,
    h = NA,
    size = NA,
    decay = NA,
    maxit = NA
  )
  
  # Initialisation des valeurs optimales pour NN seul
  optimal_nn_solo_params <- reactiveValues(
    size = NA,
    decay = NA,
    maxit = NA
  )
  
  
  
  
  # Performance analysis results
  # Bloc eventReactive optimisé
  perf_results <- eventReactive(input$run_analysis, {
    req(input$n_values)
    m_true <- parse_true_func()
    T_params <- get_t_params()
    seed <- if (is.na(input$seed)) NULL else input$seed
    
    withProgress(message = 'Running simulations', value = 0, {
      results <- list()
      n_sizes <- as.numeric(input$n_values)
      
      for (i in seq_along(n_sizes)) {
        incProgress(1 / length(n_sizes), detail = paste("n =", n_sizes[i]))
        
        # Grille pour h uniquement si sélection automatique
        h_grid <- seq(0.03, 0.4, length.out = 10)
        h_arg <- if (input$auto_h_perf) NULL else input$h_perf
        
        results[[i]] <- analyze_performance(
          n = n_sizes[i],
          N_sim = input$n_sim,
          m_true = m_true,
          T_dist = input$t_dist,
          T_params = T_params,
          h_grid = h_grid,
          seed = seed,
          h_fixed = h_arg
        )
      }
    })
    
    names(results) <- paste0("n", input$n_values)
    results
  })
  
  
  
  
  # Main plot with legend at bottom
  output$main_plot <- renderPlotly({
    results <- est_results()
    data <- results$data
    
    plot_ly() %>%
      add_trace(
        x = results$x, y = results$m_est, 
        type = 'scatter', mode = 'lines',  # Retirer 'lines+markers', ne garder que 'lines'
        name = "Estimated m̂(x)", 
        line = list(color = "#e74c3c", width = 3),
        hovertemplate = "x: %{x:.2f}<br>y: %{y:.2f}<extra></extra>"
      ) %>%
      add_trace(
        x = results$x, y = results$m_true, 
        type = 'scatter', mode = 'lines',
        name = "True m(x)", 
        line = list(color = "#3498db", dash = "dot", width = 3),
        hovertemplate = "x: %{x:.2f}<br>y: %{y:.2f}<extra></extra>"
      ) %>%
      add_markers(
        x = data$X_obs, y = data$Y_obs,
        name = "Observed Data",
        marker = list(color = "#2ecc71", size = 6, opacity = 0.5),
        hovertemplate = "X: %{x:.2f}<br>Y: %{y:.2f}<extra></extra>"
      ) %>%
      layout(
        title = list(
          text = "Nonparametric Regression with Left-Truncated Data",
          font = list(size = 20)
        ),
        xaxis = list(title = "x", gridcolor = "#e1e5e9"),
        yaxis = list(title = "m(x)", gridcolor = "#e1e5e9"),
        hoverlabel = list(bgcolor = "white", font = list(size = 14)),
        legend = list(
          orientation = "h",
          x = 0.5,
          y = -0.2,
          xanchor = "center",
          yanchor = "top"
        ),
        plot_bgcolor = "#f8f9fa",
        paper_bgcolor = "#f8f9fa",
        margin = list(t = 60, b = 100) # Increased bottom margin for legend
      ) %>%
      config(displayModeBar = TRUE)
  })
  
  output$selected_h_display <- renderPrint({
    res <- est_results()
    cat("Bande de lissage utilisée (h):", round(res$h_val, 4))
  })
  
  
  # Performance plots with legend at bottom
  output$bias_plot <- renderPlotly({
    results <- perf_results()
    if (is.null(results)) return()
    
    df <- do.call(rbind, lapply(results, function(res) {
      data.frame(
        n = res$n,
        Bias = mean(res$bias, na.rm = TRUE)
      )
    }))
    
    plot_ly(df, x = ~n, y = ~Bias, type = 'scatter', mode = 'lines+markers',
            name = "Bias", line = list(color = "red")) %>%
      layout(
        title = "Moyenne du biais selon la taille d’échantillon",
        xaxis = list(title = "Taille de l’échantillon n", type = "log"),
        yaxis = list(title = "Biais moyen"),
        plot_bgcolor = "#f8f9fa",
        legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center")
      )
  })
  
  
  
  output$variance_plot <- renderPlotly({
    results <- perf_results()
    if (is.null(results)) return()
    
    df <- do.call(rbind, lapply(results, function(res) {
      data.frame(
        n = res$n,
        Variance = mean(res$variance, na.rm = TRUE)
      )
    }))
    
    plot_ly(df, x = ~n, y = ~Variance, type = 'scatter', mode = 'lines+markers',
            name = "Variance", line = list(color = "blue")) %>%
      layout(
        title = "Moyenne de la variance selon la taille d’échantillon",
        xaxis = list(title = "Taille de l’échantillon n", type = "log"),
        yaxis = list(title = "Variance moyenne"),
        plot_bgcolor = "#f8f9fa",
        legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center")
      )
  })
  
  
  output$mse_plot <- renderPlotly({
    results <- perf_results()
    if (is.null(results)) return()
    
    df <- do.call(rbind, lapply(results, function(res) {
      data.frame(
        n = res$n,
        MSE = mean(res$mse, na.rm = TRUE)
      )
    }))
    
    plot_ly(df, x = ~n, y = ~MSE, type = 'scatter', mode = 'lines+markers',
            name = "MSE", line = list(color = "darkgreen")) %>%
      layout(
        title = "MSE selon la taille d’échantillon",
        xaxis = list(title = "Taille de l’échantillon n", type = "log"),
        yaxis = list(title = "Erreur quadratique moyenne"),
        plot_bgcolor = "#f8f9fa",
        legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center")
      )
  })
  
  
  
  # Performance metrics table
  output$metrics_table <- renderDT({
    results <- perf_results()
    if (is.null(results)) return()
    
    # Ajout de h_mean à chaque ligne
    metrics <- lapply(results, function(res) {
      data.frame(
        N = res$n,
        Avg_Bias = mean(res$bias, na.rm = TRUE),
        Avg_Variance = mean(res$variance, na.rm = TRUE),
        Avg_MSE = mean(res$mse, na.rm = TRUE),
        Mean_Bandwidth = round(res$h_mean, 4),  # 🔍 moyenne de h utilisée sur les N_sim
        Simulations = res$N_sim
      )
    })
    
    df <- do.call(rbind, metrics)
    
    datatable(
      df,
      options = list(
        pageLength = 5,
        scrollX = TRUE,
        dom = 'tip'
      ),
      rownames = FALSE,
      class = "hover stripe",
      caption = "📊 Moyennes des métriques de performance sur x et simulations"
    ) %>% formatRound(columns = 2:5, digits = 3)
  })
  
  
  # Truncation summary
  output$trunc_summary <- renderTable({
    data <- sim_data()
    total <- length(data$X_full)
    observed <- length(data$X_obs)
    
    data.frame(
      "Total Observations" = total,
      "Observed (Y ≥ T)" = observed,
      "Truncated (Y < T)" = total - observed,
      "Truncation Rate" = paste0(round(100 * (total - observed)/total, 1), "%"),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, align = "c", width = "100%")
  
  # Data table
  output$data_table <- renderDT({
    data <- sim_data()
    df <- data.frame(X = data$X_obs, Y = data$Y_obs, T = data$T_obs)
    datatable(
      df,
      options = list(
        pageLength = 5,
        scrollX = TRUE,
        dom = 'tip',
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#2c3e50', 'color': '#fff'});",
          "}")
      ),
      class = "hover stripe",
      rownames = FALSE,
      filter = "top"
    ) %>% formatRound(columns = 1:3, digits = 3)
  })
  
  # Bandwidth diagnostic plot
  output$bw_diagnostic_plot <- renderPlotly({
    req(input$h_range)
    data <- sim_data()
    m_true <- parse_true_func()
    x_grid <- seq(0, 1, length.out = 100)
    h_values <- seq(input$h_range[1], input$h_range[2], length.out = 20)
    
    mse_values <- sapply(h_values, function(h) {
      m_est <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h)
      if (any(is.na(m_est))) return(NA)
      mean((m_est - m_true(x_grid))^2, na.rm = TRUE)
    })
    
    plot_ly(x = h_values, y = mse_values, type = 'scatter', mode = 'lines+markers',
            line = list(color = '#9b59b6', width = 3),
            marker = list(color = '#9b59b6', size = 8),
            name = "MSE vs Bandwidth") %>%
      add_trace(x = c(input$h, input$h), y = range(mse_values, na.rm = TRUE),
                type = 'scatter', mode = 'lines',
                line = list(color = '#e74c3c', dash = 'dash'),
                name = "Current h") %>%
      layout(
        xaxis = list(title = "Bandwidth h", gridcolor = "#e1e5e9"),
        yaxis = list(title = "Mean Squared Error", gridcolor = "#e1e5e9"),
        title = "Bandwidth Selection Diagnostic",
        plot_bgcolor = "#f8f9fa",
        hovermode = "x unified",
        legend = list(
          orientation = "h",
          x = 0.5,
          y = -0.3,
          xanchor = "center",
          yanchor = "top"
        ),
        margin = list(b = 100)
      )
  })
  
  
  
  # Bandwidth study results
  h_study_results <- eventReactive(input$run_h_study, {
    req(input$h_range_study, input$h_steps)
    m_true <- parse_true_func()
    T_params <- get_t_params()
    seed <- if (is.na(input$seed)) NULL else input$seed
    
    h_values <- seq(input$h_range_study[1], input$h_range_study[2], 
                    length.out = input$h_steps)
    
    withProgress(message = 'Running bandwidth study', value = 0, {
      data <- simulate_data(input$h_study_n, m_true, input$t_dist, T_params, seed)
      x_grid <- seq(0, 1, length.out = 200)
      
      results <- list()
      for (i in seq_along(h_values)) {
        incProgress(1/length(h_values), detail = paste("h =", round(h_values[i], 3)))
        results[[i]] <- list(
          h = h_values[i],
          m_est = nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h_values[i])
        )
      }
    })
    
    list(
      x_grid = x_grid,
      m_true = m_true(x_grid),
      data = data,
      results = results
    )
  })
  
  output$h_study_plot <- renderPlotly({
    results <- h_study_results()
    if (is.null(results)) return()
    
    # Préparation des frames
    frames <- lapply(results$results, function(res) {
      list(
        name = sprintf("h%.3f", res$h),
        data = list(
          list(x = results$x_grid, y = results$m_true,
               type = 'scatter', mode = 'lines',
               line = list(dash = "dash", color = "blue"),
               name = "True m(x)"),
          list(x = results$x_grid, y = res$m_est,
               type = 'scatter', mode = 'lines',
               line = list(color = "red"),
               name = "Estimated m̂(x)")
        )
      )
    })
    
    h_values <- sapply(results$results, function(x) x$h)
    
    # Graphique initial avec slider seulement
    fig <- plot_ly() %>%
      add_trace(
        x = results$x_grid, y = results$m_true,
        type = 'scatter', mode = 'lines',
        line = list(dash = "dash", color = "blue"),
        name = "True m(x)"
      ) %>%
      add_trace(
        x = results$x_grid, y = results$results[[1]]$m_est,
        type = 'scatter', mode = 'lines',
        line = list(color = "red"),
        name = "Estimated m̂(x)"
      ) %>%
      layout(
        title = "Effect of Bandwidth Parameter h",
        xaxis = list(title = "x"),
        yaxis = list(title = "m(x)"),
        sliders = list(
          list(
            active = 0,
            currentvalue = list(prefix = "h = "),
            steps = lapply(seq_along(h_values), function(i) {
              list(
                label = sprintf("%.3f", h_values[i]),
                method = "animate",
                args = list(
                  list(sprintf("h%.3f", h_values[i])),
                  list(mode = "immediate", frame = list(duration = 0), transition = list(duration = 0))
                )
              )
            })
          )
        ),
        legend = list(
          orientation = "h",
          x = 0.5,
          y = -0.2,
          xanchor = "center",
          yanchor = "top"
        ),
        margin = list(b = 100)
      )
    
    fig$x$frames <- frames
    fig
  })
  
  
  
  
  
  
  
  
  
  # Multiple bandwidth plots
  output$h_multiple_plots <- renderPlot({
    results <- h_study_results()
    if (is.null(results)) return()
    
    # Select a subset of bandwidths to display (max 4)
    h_values <- sapply(results$results, function(x) x$h)
    idx <- seq(1, length(h_values), length.out = min(4, length(h_values)))
    selected_results <- results$results[idx]
    
    plots <- lapply(selected_results, function(res) {
      df <- data.frame(
        x = results$x_grid,
        m_true = results$m_true,
        m_est = res$m_est
      )
      
      ggplot(df, aes(x = x)) +
        geom_line(aes(y = m_true, color = "True m(x)"), linetype = "dashed", size = 1) +
        geom_line(aes(y = m_est, color = "Estimated m̂(x)"), size = 1) +
        scale_color_manual(values = c("True m(x)" = "blue", "Estimated m̂(x)" = "red")) +
        labs(title = paste("h =", round(res$h, 3)), y = "m(x)", x = "x", color = "") +
        theme_minimal() +
        theme(legend.position = "bottom")
    })
    
    # Arrange plots in a grid
    grid.arrange(grobs = plots, ncol = 2)
  })
  
  # Dans le server
  output$download_h_multiple <- downloadHandler(
    filename = function() {
      paste("bandwidth_comparison", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      results <- h_study_results()
      if (is.null(results)) return()
      
      h_values <- sapply(results$results, function(x) x$h)
      idx <- seq(1, length(h_values), length.out = min(4, length(h_values)))
      selected_results <- results$results[idx]
      
      plots <- lapply(selected_results, function(res) {
        df <- data.frame(
          x = results$x_grid,
          m_true = results$m_true,
          m_est = res$m_est
        )
        
        ggplot(df, aes(x = x)) +
          geom_line(aes(y = m_true, color = "True m(x)"), linetype = "dashed", size = 1) +
          geom_line(aes(y = m_est, color = "Estimated m̂(x)"), size = 1) +
          scale_color_manual(values = c("True m(x)" = "blue", "Estimated m̂(x)" = "red")) +
          labs(title = paste("h =", round(res$h, 3)), y = "m(x)", x = "x", color = "") +
          theme_minimal() +
          theme(legend.position = "bottom")
      })
      
      combined_plot <- arrangeGrob(grobs = plots, ncol = 2)
      ggsave(file, combined_plot, width = 10, height = 8)
    }
  )
  
  # MISE simulation results
  mise_results <- eventReactive(input$run_mise_simulation, {
    req(input$mise_n_values)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    seed <- if (is.na(input$seed)) NULL else input$seed
    
    n_vals <- seq(input$mise_n_values[1], input$mise_n_values[2], length.out = 3)
    h_grid <- seq(0.03, 0.4, length.out = 10)
    
    withProgress(message = 'Running MISE simulation', value = 0, {
      results <- lapply(n_vals, function(n_val) {
        incProgress(1 / length(n_vals), detail = paste("n =", n_val))
        
        res <- compute_mise(
          N = n_val,
          m_true = m_true,
          T_dist = input$t_dist,
          T_params = T_params,
          N_sim = input$mise_n_sim,
          h_grid = h_grid,
          seed = seed
        )
        
        data.frame(n = res$n, mise = res$mise, h_mean = res$h_mean)
      })
      
      do.call(rbind, results)
    })
  })
  
  
  
  # MISE line plot modifié
  output$mise_line_plot <- renderPlotly({
    results <- mise_results()
    if (is.null(results)) return()
    
    plot_ly(results, x = ~n, y = ~mise, type = 'scatter', mode = 'lines+markers',
            text = ~paste("h̄ =", round(h_mean, 3), "<br>MISE:", round(mise, 4))) %>%
      layout(
        title = "MISE en fonction de n (avec hₙ optimal)",
        xaxis = list(title = "Taille de l'échantillon (n)"),
        yaxis = list(title = "MISE"),
        plot_bgcolor = "#f8f9fa"
      )
  })
  
  
  output$mise_table <- renderDT({
    results <- mise_results()
    if (is.null(results)) return()
    print("🧪 Colnames de mise_results() :")
    print(colnames(results))
    
    print("🧪 Aperçu des données :")
    print(head(results, 3))
    
    datatable(
      results,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    ) %>% formatRound(columns = intersect(c("mise", "h_mean"), colnames(results)), digits = 3)
  })
  
  
  
  
  
  
  
  
  
  # Comparison of estimators
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
      if (input$auto_h_comp) {
        h_cv <- select_h_cv(data$X_obs, data$Y_obs, data$T_obs, seq(0.03, 0.4, length.out = 10))
      } else {
        h_cv <- input$comp_h
      }
      nw_est <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h_cv)
      
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
      
      
      
      incProgress(0.2, detail = "Finalizing results")
      
      list(
        x_grid = x_grid,
        m_true = m_vals,
        results = results,
        data = data
      )
    })
  })
  
  # Comparison plot
  output$comparison_plot <- renderPlot({
    results <- comparison_results()
    if (is.null(results)) return()
    
    # Create data frames for each estimator
    df_nw <- data.frame(x = results$x_grid, 
                        True = results$m_true, 
                        Estimate = results$results$nw$estimate,
                        Method = "Nadaraya-Watson")
    
    df_loess <- data.frame(x = results$x_grid, 
                           True = results$m_true, 
                           Estimate = results$results$loess$estimate,
                           Method = "LOESS")
    
    df_poly <- data.frame(x = results$x_grid, 
                          True = results$m_true, 
                          Estimate = results$results$poly$estimate,
                          Method = "Local Polynomial")
    
    df_spline <- data.frame(x = results$x_grid, 
                            True = results$m_true, 
                            Estimate = results$results$spline$estimate,
                            Method = "Splines")
    
    # Combine all data
    df_all <- rbind(df_nw, df_loess, df_poly, df_spline)
    
    # Create plots
    p <- ggplot(df_all, aes(x = x)) +
      geom_line(aes(y = True, color = "True Function"), size = 1, linetype = "dashed") +
      geom_line(aes(y = Estimate, color = Method), size = 1) +
      facet_wrap(~Method, ncol = 2) +
      scale_color_manual(values = c("True Function" = "black",
                                    "Nadaraya-Watson" = "red",
                                    "LOESS" = "blue",
                                    "Local Polynomial" = "darkgreen",
                                    "Splines" = "purple")) +
      labs(title = "Comparison of Nonparametric Estimators with Truncation",
           x = "x", y = "m(x)", color = "") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    p
  })
  
  output$download_comparison_plot <- downloadHandler(
    filename = function() {
      paste("comparison_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      results <- comparison_results()
      if (is.null(results)) return()
      
      df_nw <- data.frame(x = results$x_grid, True = results$m_true, 
                          Estimate = results$results$nw$estimate, Method = "Nadaraya-Watson")
      df_loess <- data.frame(x = results$x_grid, True = results$m_true, 
                             Estimate = results$results$loess$estimate, Method = "LOESS")
      df_poly <- data.frame(x = results$x_grid, True = results$m_true, 
                            Estimate = results$results$poly$estimate, Method = "Local Polynomial")
      df_spline <- data.frame(x = results$x_grid, True = results$m_true, 
                              Estimate = results$results$spline$estimate, Method = "Splines")
      df_all <- rbind(df_nw, df_loess, df_poly, df_spline)
      
      p <- ggplot(df_all, aes(x = x)) +
        geom_line(aes(y = True, color = "True Function"), size = 1, linetype = "dashed") +
        geom_line(aes(y = Estimate, color = Method), size = 1) +
        facet_wrap(~Method, ncol = 2) +
        scale_color_manual(values = c("True Function" = "black",
                                      "Nadaraya-Watson" = "red",
                                      "LOESS" = "blue",
                                      "Local Polynomial" = "darkgreen",
                                      "Splines" = "purple")) +
        labs(title = "Comparison of Nonparametric Estimators with Truncation",
             x = "x", y = "m(x)", color = "") +
        theme_minimal() +
        theme(legend.position = "bottom",
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA))
      
      ggsave(file, p, width = 10, height = 8)
    }
  )
  
  
  # MSE comparison plot
  output$mse_comparison_plot <- renderPlotly({
    results <- comparison_results()
    if (is.null(results)) return()
    
    mse_data <- data.frame(
      Method = c("Nadaraya-Watson", "LOESS", "Local Polynomial", "Splines"),
      MSE = c(results$results$nw$mse,
              results$results$loess$mse,
              results$results$poly$mse,
              results$results$spline$mse)
    )
    
    plot_ly(mse_data, x = ~Method, y = ~MSE, type = 'bar',
            color = ~Method,
            colors = c("Nadaraya-Watson" = "red",
                       "LOESS" = "blue",
                       "Local Polynomial" = "darkgreen",
                       "Splines" = "purple")) %>%
      layout(title = "Mean Squared Error Comparison",
             xaxis = list(title = "Estimator"),
             yaxis = list(title = "MSE"),
             showlegend = FALSE)
  })
  
  
  
  
  
  
  # Uniform convergence analysis
  unif_results <- eventReactive(input$run_unif_convergence, {
    req(input$unif_n_values)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    seed <- if (is.na(input$seed)) NULL else input$seed
    
    n_vec <- as.numeric(input$unif_n_values)
    x_grid <- seq(0.05, 0.95, length.out = 100)
    true_vals <- m_true(x_grid)
    h_grid <- seq(0.03, 0.4, length.out = 10)  # grille utilisée pour CV
    
    withProgress(message = 'Analyzing uniform convergence', value = 0, {
      sup_errors <- matrix(NA, nrow = length(n_vec), ncol = input$unif_n_sim)
      
      for (i in seq_along(n_vec)) {
        n <- n_vec[i]
        
        for (s in 1:input$unif_n_sim) {
          incProgress(1 / (length(n_vec) * input$unif_n_sim),
                      detail = paste("n =", n, "sim", s, "of", input$unif_n_sim))
          
          if (!is.null(seed)) set.seed(seed + i * input$unif_n_sim + s)
          data <- simulate_data(n, m_true, input$t_dist, T_params)
          
          h_val <- if (input$auto_h_unif) {
            select_h_cv(data$X_obs, data$Y_obs, data$T_obs, h_grid)
          } else {
            input$unif_h
          }
          
          m_hat <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h_val)
          sup_errors[i, s] <- max(abs(m_hat - true_vals), na.rm = TRUE)
        }
      }
      
      mean_sup_errors <- rowMeans(sup_errors, na.rm = TRUE)
      
      list(
        n_vec = n_vec,
        mean_errors = mean_sup_errors,
        error_matrix = sup_errors,
        n_sim = input$unif_n_sim
      )
    })
  })
  
  
  
  # Uniform convergence plot
  output$unif_convergence_plot <- renderPlotly({
    results <- unif_results()
    if (is.null(results)) return()
    
    df <- data.frame(
      n = results$n_vec,
      mean_error = results$mean_errors
    )
    
    plot_ly(df, x = ~n, y = ~mean_error, type = 'scatter', mode = 'lines+markers',
            line = list(color = 'darkgreen', width = 2),
            marker = list(color = 'darkgreen', size = 8),
            text = ~paste("n =", n, "<br>Mean error:", round(mean_error, 4)),
            hoverinfo = 'text') %>%
      layout(
        title = "Uniform Almost Sure Convergence",
        xaxis = list(title = "Sample size n", type = "log"),
        yaxis = list(title = "Mean maximum absolute error"),
        showlegend = FALSE
      )
  })
  
  # Error matrix table
  output$unif_error_table <- renderDT({
    results <- unif_results()
    if (is.null(results)) return()
    
    # Create a summary data frame
    df <- data.frame(
      SampleSize = results$n_vec,
      MeanError = results$mean_errors,
      MinError = apply(results$error_matrix, 1, min, na.rm = TRUE),
      MaxError = apply(results$error_matrix, 1, max, na.rm = TRUE),
      StdDev = apply(results$error_matrix, 1, sd, na.rm = TRUE)
    )
    
    datatable(
      df,
      options = list(
        pageLength = 5,
        scrollX = TRUE,
        dom = 'tip'
      ),
      rownames = FALSE,
      class = "hover stripe",
      caption = paste("Uniform error statistics from", results$n_sim, "simulations"),
      colnames = c('Sample Size', 'Mean Error', 'Minimum Error', 'Maximum Error', 'Std. Deviation')
    ) %>% formatRound(columns = 2:5, digits = 3)
  })
  
  
  
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
  
  
  
  # Neural network results
  # Neural network results
  nn_results <- eventReactive(input$run_mix_plot, {
    req(input$nn_sample_size)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    seed <- input$nn_seed
    if (!is.null(seed)) set.seed(seed)
    
    data <- simulate_data(
      n = input$nn_sample_size,
      m_true = m_true,
      t_dist = input$t_dist,
      t_params = T_params,
      seed = seed
    )
    
    x_grid <- seq(0, 1, length.out = 100)
    
    # Sélection de h (auto ou fixe)
    h_val <- if (isTRUE(input$auto_h_mix)) {
      select_h_cv(data$X_obs, data$Y_obs, data$T_obs, h_grid = seq(0.03, 0.4, length.out = 10))
    } else {
      input$mix_h
    }
    
    nn_result <- train_nn_estimator(
      X_obs = data$X_obs,
      Y_obs = data$Y_obs,
      T_obs = data$T_obs,
      size = input$nn_size,
      decay = input$nn_decay,
      maxit = input$nn_maxit,
      x_grid = x_grid
    )
    
    list(
      result = nn_result,
      x_grid = x_grid,
      m_true = m_true(x_grid),
      data = data,
      h_val = h_val
    )
  })
  
  
  output$nn_mix_comparison_plot <- renderPlotly({
    req(input$run_mix_plot)
    results <- nn_results()
    if (is.null(results)) return()
    
    lambda <- input$mix_lambda
    h_val <- results$h_val
    
    x <- results$x_grid
    m_true <- results$m_true
    nw <- nw_estimator_truncated(x, results$data$X_obs, results$data$Y_obs, results$data$T_obs, h_val)
    nn <- results$result$predictions
    mix <- lambda * nw + (1 - lambda) * nn
    
    df_all <- rbind(
      data.frame(x = x, y = m_true, label = "True", group = "NW vs True"),
      data.frame(x = x, y = nw, label = "NW", group = "NW vs True"),
      
      data.frame(x = x, y = m_true, label = "True", group = "NN vs True"),
      data.frame(x = x, y = nn, label = "NN", group = "NN vs True"),
      
      data.frame(x = x, y = m_true, label = "True", group = "Mixture vs True"),
      data.frame(x = x, y = mix, label = "Mixture", group = "Mixture vs True")
    )
    
    p <- ggplot(df_all, aes(x = x, y = y, color = label, linetype = label)) +
      geom_line(size = 1) +
      facet_wrap(~group, ncol = 1, scales = "free_y") +
      labs(
        title = "Comparaison par méthode avec la fonction vraie",
        x = "x", y = "m(x)", color = "Courbe", linetype = "Style"
      ) +
      scale_color_manual(values = c("True" = "blue", "NW" = "green", "NN" = "red", "Mixture" = "purple")) +
      scale_linetype_manual(values = c("True" = "dashed", "NW" = "solid", "NN" = "solid", "Mixture" = "solid")) +
      theme_minimal()
    
    ggplotly(p)
  })
  
  
  
  
  output$y_density_plot <- renderPlotly({
    data <- sim_data()
    
    df <- data.frame(
      Y = c(data$Y_full, data$Y_obs),
      Source = factor(rep(c("Avant troncature", "Après troncature"),
                          c(length(data$Y_full), length(data$Y_obs))))
    )
    
    p <- ggplot(df, aes(x = Y, fill = Source, color = Source)) +
      geom_density(alpha = 0.4, size = 1.2) +
      labs(title = "Densité de Y avant et après troncature",
           x = "Valeur de Y", y = "Densité") +
      scale_fill_manual(values = c("Avant troncature" = "#3498db", "Après troncature" = "#e74c3c")) +
      scale_color_manual(values = c("Avant troncature" = "#3498db", "Après troncature" = "#e74c3c")) +
      theme_minimal(base_size = 15) +
      theme(legend.position = "bottom")
    
    ggplotly(p, tooltip = c("x", "y", "fill"))
  })
  
  
  
  
  output$mse_vs_h_plot <- renderPlotly({
    results <- h_study_results()
    if (is.null(results)) return()
    
    mse_vals <- sapply(results$results, function(res) {
      mean((res$m_est - results$m_true)^2, na.rm = TRUE)
    })
    
    df <- data.frame(h = sapply(results$results, `[[`, "h"), mse = mse_vals)
    
    plot_ly(df, x = ~h, y = ~mse, type = "scatter", mode = "lines+markers",
            line = list(color = "#2980b9"), marker = list(size = 6)) %>%
      layout(
        title = "MSE en fonction du paramètre de lissage h",
        xaxis = list(title = "h"),
        yaxis = list(title = "Mean Squared Error"),
        plot_bgcolor = "#f8f9fa"
      )
  })
  
  
  output$bandwidth_table <- DT::renderDT({
    results <- h_study_results()
    if (is.null(results)) return()
    
    df <- data.frame(
      h = round(sapply(results$results, `[[`, "h"), 3),
      MSE = round(sapply(results$results, function(res) mean((res$m_est - results$m_true)^2, na.rm = TRUE)), 3)
    )
    
    datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  
  output$h_best_summary <- renderPrint({
    results <- h_study_results()
    if (is.null(results)) return()
    
    mse_vals <- sapply(results$results, function(res) mean((res$m_est - results$m_true)^2, na.rm = TRUE))
    h_vals <- sapply(results$results, `[[`, "h")
    h_opt <- h_vals[which.min(mse_vals)]
    
    cat("✅ Le h optimal (minimisant la MSE) est :", round(h_opt, 4), "\n")
  })
  
  
  
  output$mise_best_summary <- renderPrint({
    results <- mise_results()
    if (is.null(results)) return()
    
    best <- results[which.min(results$mise), ]
    cat("✅ Meilleure MISE obtenue pour :\n")
    cat("- n =", best$n, "\n")
    cat("- h̄ (moyenne des hₙ optimaux) =", round(best$h_mean,3), "\n")
    cat("- MISE =", round(best$mise, 3), "\n")
  })
  
  
  
  
  
  
  
  
  output$loglog_unif_plot <- renderPlotly({
    results <- unif_results()
    if (is.null(results)) return()
    
    log_n <- log(results$n_vec)
    log_e <- log(results$mean_errors)
    model <- lm(log_e ~ log_n)
    
    df <- data.frame(log_n = log_n, log_e = log_e)
    
    plot_ly(df, x = ~log_n, y = ~log_e, type = "scatter", mode = "markers",
            marker = list(size = 8, color = "#2980b9")) %>%
      add_lines(x = ~log_n, y = fitted(model), line = list(color = "red"), name = "Régression") %>%
      layout(
        title = "Log-Log Régression : vitesse de convergence",
        xaxis = list(title = "log(n)"),
        yaxis = list(title = "log(sup error)"),
        annotations = list(
          list(
            x = 0.1, y = 0.95, xref = "paper", yref = "paper",
            text = paste0("Pente estimée ≈ ", round(coef(model)[2], 3)),
            showarrow = FALSE, font = list(size = 12, color = "black")
          )
        )
      )
  })
  
  
  
  output$best_estimator_summary <- renderPrint({
    res <- efficiency_results()
    if (is.null(res)) return()
    best <- names(which.min(res$mse))
    cat("✅ Best estimator (lowest MSE):", best)
  })
  
  
  # EventReactive pour lancer l'étude
  lambda_study_results <- eventReactive(input$run_lambda_study, {
    req(input$lambda_n, input$lambda_nn_size, input$lambda_nn_decay, input$lambda_nn_maxit, input$lambda_nsim)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    seed <- input$lambda_seed
    if (!is.null(seed)) set.seed(seed)
    
    lambda_vals <- seq(0, 1, by = 0.1)
    x_grid <- seq(0, 1, length.out = 20)
    h_grid <- seq(0.03, 0.4, length.out = 10)
    
    mse_lambda <- numeric(length(lambda_vals))
    
    withProgress(message = "Simulation selon λ...", value = 0, {
      for (i in seq_along(lambda_vals)) {
        lambda <- lambda_vals[i]
        mse_list <- numeric(input$lambda_nsim)
        
        for (b in 1:input$lambda_nsim) {
          if (!is.null(seed)) set.seed(seed + b)
          data <- simulate_data(input$lambda_n, m_true, input$t_dist, T_params)
          
          # h_n automatique ou fixe
          h_val <- if (isTRUE(input$auto_h_lambda)) {
            select_h_cv(data$X_obs, data$Y_obs, data$T_obs, h_grid = h_grid)
          } else {
            input$lambda_h
          }
          
          nw <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h_val)
          
          nn_result <- train_nn_estimator(
            X_obs = data$X_obs,
            Y_obs = data$Y_obs,
            T_obs = data$T_obs,
            size = input$lambda_nn_size,
            decay = input$lambda_nn_decay,
            maxit = input$lambda_nn_maxit,
            x_grid = x_grid
          )
          
          mix <- lambda * nw + (1 - lambda) * nn_result$predictions
          mse_list[b] <- mean((mix - m_true(x_grid))^2, na.rm = TRUE)
        }
        
        mse_lambda[i] <- mean(mse_list, na.rm = TRUE)
        incProgress(1 / length(lambda_vals), detail = paste("λ =", lambda))
      }
    })
    
    # Filtrer pour ne conserver que les lambda strictement positifs
    valid_indices <- which(lambda_vals > 0)
    if (length(valid_indices) == 0) {
      warning("Aucune valeur de lambda strictement positive trouvée. Vérifie la définition de lambda_vals.")
      return(list(lambda_vals = lambda_vals, mse = mse_lambda, best_lambda = NA))
    }
    
    lambda_vals <- lambda_vals[valid_indices]
    mse_lambda <- mse_lambda[valid_indices]
    
    best_lambda <- lambda_vals[which.min(mse_lambda)]
    optimal_params$lambda <- best_lambda
    # Si h est automatique, stockons aussi la dernière valeur de h utilisée
    if (isTRUE(input$auto_h_lambda)) {
      optimal_params$h <- h_val
    } else {
      optimal_params$h <- input$lambda_h
    }
    
    list(
      lambda_vals = lambda_vals,
      mse = mse_lambda,
      best_lambda = best_lambda
    )
  })
  
  
  
  
  output$lambda_mse_plot <- renderPlotly({
    res <- lambda_study_results()
    if (is.null(res)) return()
    
    df <- data.frame(lambda = res$lambda_vals, MSE = res$mse)
    
    plot_ly(df, x = ~lambda, y = ~MSE, type = 'scatter', mode = 'lines+markers',
            name = "MSE(λ)", line = list(color = "blue")) %>%
      add_trace(x = res$best_lambda, y = min(df$MSE),
                type = "scatter", mode = "markers+text",
                text = paste0("λ* = ", round(res$best_lambda, 2)),
                textposition = "top right",
                marker = list(color = "red", size = 10),
                name = "λ optimal") %>%
      layout(
        title = "Évolution de l'erreur quadratique moyenne selon λ",
        xaxis = list(title = "λ (Poids pour NW)"),
        yaxis = list(title = "MSE moyen"),
        legend = list(
          orientation = "h", x = 0.5, y = -0.25,
          xanchor = "center", font = list(size = 13)
        ),
        plot_bgcolor = "#f9f9f9"
      )
  })
  
  output$best_lambda_display <- renderPrint({
    res <- lambda_study_results()
    if (is.null(res)) return()
    cat("✅ λ optimal (minimisant le MSE) :", round(res$best_lambda, 3))
  })
  
  
  nn_param_study_results <- eventReactive(input$run_nn_param_study, {
    req(input$nn_study_nsim, input$nn_study_lambda)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    x_grid <- seq(0, 1, length.out = 50)
    h_grid <- seq(0.03, 0.4, length.out = 10)
    
    param_vals <- switch(input$nn_param_to_study,
                         "size" = c(1, 5, 10, 20, 30, 50),
                         "decay" = c(0, 0.001, 0.01, 0.05, 0.1),
                         "maxit" = c(100, 500, 1000, 3000)
    )
    
    grid_results <- list()
    
    withProgress(message = "⏳ NN parameters", value = 0, {
      for (i in seq_along(param_vals)) {
        pval <- param_vals[i]
        mse_vec <- numeric(input$nn_study_nsim)
        
        for (b in 1:input$nn_study_nsim) {
          data <- simulate_data(input$nn_study_n, m_true, input$t_dist, T_params)
          
          h_val <- if (isTRUE(input$auto_h_nn_param)) {
            select_h_cv(data$X_obs, data$Y_obs, data$T_obs, h_grid = h_grid)
          } else {
            input$nn_study_h
          }
          
          nw <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h_val)
          nn <- train_nn_estimator(
            X_obs = data$X_obs, Y_obs = data$Y_obs, T_obs = data$T_obs,
            size = if (input$nn_param_to_study == "size") pval else input$nn_fixed_size,
            decay = if (input$nn_param_to_study == "decay") pval else input$nn_fixed_decay,
            maxit = if (input$nn_param_to_study == "maxit") pval else input$nn_fixed_maxit,
            x_grid = x_grid
          )$predictions
          
          mix <- input$nn_study_lambda * nw + (1 - input$nn_study_lambda) * nn
          mse_vec[b] <- mean((mix - m_true(x_grid))^2, na.rm = TRUE)
        }
        
        grid_results[[i]] <- data.frame(
          Param = pval,
          MSE = mean(mse_vec)
        )
        
        incProgress(1 / length(param_vals), detail = paste("Param =", pval))
      }
    })
    
    df_all <- do.call(rbind, grid_results)
    
    # 🔥 Filtrer pour ne conserver que les paramètres strictement positifs
    df_all <- df_all[df_all$Param > 0, ]
    
    if (nrow(df_all) == 0) {
      warning("Aucun paramètre strictement positif trouvé. Vérifie la définition de param_vals.")
      return(list(df = data.frame(), best_param = NA, best_mse = NA))
    }
    
    best_row <- df_all[which.min(df_all$MSE), ]
    best_param <- best_row$Param
    optimal_params[[input$nn_param_to_study]] <- best_param
    
    list(
      df = df_all,
      best_param = best_param,
      best_mse = best_row$MSE
    )
  })
  
  
  
  
  
  
  
  
  
  
  output$nn_param_study_plot <- renderPlotly({
    res <- nn_param_study_results()
    if (is.null(res)) return()
    df <- res$df
    
    plot_ly(df, x = ~Param) %>%
      add_trace(y = ~MSE, type = 'scatter', mode = 'lines+markers',
                name = "MSE", line = list(color = "darkgreen")) %>%
      add_trace(
        x = res$best_param, y = res$best_mse, mode = "markers+text",
        text = paste0("Optimal = ", res$best_param),
        marker = list(color = "red", size = 10),
        name = "Optimal"
      ) %>%
      layout(
        title = paste("MSE moyen selon", input$nn_param_to_study),
        xaxis = list(title = input$nn_param_to_study),
        yaxis = list(title = "MSE moyen"),
        legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center")
      )
  })
  
  
  output$best_nn_param_display <- renderPrint({
    res <- nn_param_study_results()
    if (is.null(res)) return()
    paste("✓ Meilleur paramètre", input$nn_param_to_study, ":", res$best_param)
  })
  
  output$summary_optimal_params <- DT::renderDataTable({
    df <- data.frame(
      `Paramètre` = c("λ (lambda)", "h (bande de lissage)", "size (neurones)", "decay (weight decay)", "maxit (itérations)"),
      `Valeur optimale` = c(
        if (!is.null(optimal_params$lambda) && !is.na(optimal_params$lambda)) round(optimal_params$lambda, 3) else "Non défini",
        if (!is.null(optimal_params$h) && !is.na(optimal_params$h)) round(optimal_params$h, 3) else "Non défini",
        if (!is.null(optimal_params$size) && !is.na(optimal_params$size)) optimal_params$size else "Non défini",
        if (!is.null(optimal_params$decay) && !is.na(optimal_params$decay)) optimal_params$decay else "Non défini",
        if (!is.null(optimal_params$maxit) && !is.na(optimal_params$maxit)) optimal_params$maxit else "Non défini"
      )
    )
    DT::datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 5))
  })
  
  
  
  
  nn_solo_param_results <- eventReactive(input$run_nn_solo_study, {
    req(input$nn_solo_nsim, input$nn_solo_n)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    x_grid <- seq(0, 1, length.out = 100)
    
    param_vals <- switch(input$nn_solo_param,
                         "size" = c(1, 5, 10, 20, 30, 50),
                         "decay" = c(0, 0.001, 0.01, 0.05, 0.1),
                         "maxit" = c(100, 500, 1000, 3000))
    
    grid_results <- list()
    
    withProgress(message = "Étude NN seul en cours", value = 0, {
      for (i in seq_along(param_vals)) {
        pval <- param_vals[i]
        mse_vec <- numeric(input$nn_solo_nsim)
        
        for (b in 1:input$nn_solo_nsim) {
          data <- simulate_data(input$nn_solo_n, m_true, input$t_dist, T_params)
          nn <- train_nn_estimator(
            X_obs = data$X_obs, Y_obs = data$Y_obs, T_obs = data$T_obs,
            size = if (input$nn_solo_param == "size") pval else input$nn_solo_fixed_size,
            decay = if (input$nn_solo_param == "decay") pval else input$nn_solo_fixed_decay,
            maxit = if (input$nn_solo_param == "maxit") pval else input$nn_solo_fixed_maxit,
            x_grid = x_grid
          )$predictions
          
          mse_vec[b] <- mean((nn - m_true(x_grid))^2, na.rm = TRUE)
        }
        
        grid_results[[i]] <- data.frame(Param = pval, MSE = mean(mse_vec))
        incProgress(1 / length(param_vals), detail = paste(input$nn_solo_param, ":", pval))
      }
    })
    
    df_all <- do.call(rbind, grid_results)
    
    # 🔥 Filtrer pour ne conserver que les paramètres strictement positifs
    df_all <- df_all[df_all$Param > 0, ]
    
    if (nrow(df_all) == 0) {
      warning("Aucun paramètre strictement positif trouvé. Vérifie la définition de param_vals.")
      return(list(df = data.frame(), best_param = NA, best_mse = NA))
    }
    
    best_row <- df_all[which.min(df_all$MSE), ]
    best_param <- best_row$Param
    
    optimal_nn_solo_params[[input$nn_solo_param]] <- best_param
    
    list(
      df = df_all,
      best_param = best_param,
      best_mse = best_row$MSE
    )
  })
  
  
  output$nn_solo_param_plot <- renderPlotly({
    res <- nn_solo_param_results()
    if (is.null(res)) return()
    df <- res$df
    
    plot_ly(df, x = ~Param) %>%
      add_trace(y = ~MSE, type = 'scatter', mode = 'lines+markers', name = "MSE", line = list(color = "darkgreen")) %>%
      add_trace(
        x = res$best_param, y = res$best_mse, mode = "markers+text",
        text = paste0("Optimal = ", res$best_param),
        marker = list(color = "red", size = 10), name = "Optimal"
      ) %>%
      layout(
        title = paste("MSE selon", input$nn_solo_param),
        xaxis = list(title = input$nn_solo_param),
        yaxis = list(title = "MSE moyen"),
        legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center")
      )
  })
  
  output$best_nn_solo_param_display <- renderPrint({
    res <- nn_solo_param_results()
    if (is.null(res)) return()
    paste("✓ Meilleur paramètre", input$nn_solo_param, ":", res$best_param)
  })
  
  output$summary_nn_solo_params <- DT::renderDataTable({
    df <- data.frame(
      `Paramètre` = c("size (neurones)", "decay (weight decay)", "maxit (itérations)"),
      `Valeur optimale` = c(
        ifelse(is.na(optimal_nn_solo_params$size), "Non défini", optimal_nn_solo_params$size),
        ifelse(is.na(optimal_nn_solo_params$decay), "Non défini", optimal_nn_solo_params$decay),
        ifelse(is.na(optimal_nn_solo_params$maxit), "Non défini", optimal_nn_solo_params$maxit)
      )
    )
    DT::datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 5))
  })
  
  
  final_comparison_results <- eventReactive(input$run_final_comparison, {
    req(input$final_n, input$final_nsim)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    x_grid <- seq(0, 1, length.out = 100)
    
    mse_nw <- mse_nn <- mse_mix <- numeric(input$final_nsim)
    
    withProgress(message = "Comparaison en cours", value = 0, {
      for (b in 1:input$final_nsim) {
        data <- simulate_data(input$final_n, m_true, input$t_dist, T_params)
        
        # NW seul
        nw_est <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, input$final_h_nw)
        mse_nw[b] <- mean((nw_est - m_true(x_grid))^2, na.rm = TRUE)
        
        # NN seul
        nn_pred <- train_nn_estimator(
          X_obs = data$X_obs, Y_obs = data$Y_obs, T_obs = data$T_obs,
          size = input$final_size_nn, decay = input$final_decay_nn,
          maxit = input$final_maxit_nn, x_grid = x_grid
        )$predictions
        mse_nn[b] <- mean((nn_pred - m_true(x_grid))^2, na.rm = TRUE)
        
        # Mixture
        nn_mix_pred <- train_nn_estimator(
          X_obs = data$X_obs, Y_obs = data$Y_obs, T_obs = data$T_obs,
          size = input$final_size_nn_mix, decay = input$final_decay_nn_mix,
          maxit = input$final_maxit_nn_mix, x_grid = x_grid
        )$predictions
        nw_mix <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, input$final_h_mix)
        mix_pred <- input$final_lambda_mix * nw_mix + (1 - input$final_lambda_mix) * nn_mix_pred
        mse_mix[b] <- mean((mix_pred - m_true(x_grid))^2, na.rm = TRUE)
        
        incProgress(1 / input$final_nsim)
      }
    })
    
    data.frame(
      Estimateur = c("Nadaraya-Watson", "Neural Network", "Mixture NW+NN"),
      MSE = round(c(mean(mse_nw), mean(mse_nn), mean(mse_mix)), 5)
    )
  })
  
  
  output$final_comparison_plot <- renderPlotly({
    df <- final_comparison_results()
    if (is.null(df)) return()
    
    plot_ly(df, x = ~Estimateur, y = ~MSE, type = 'bar', name = "MSE",
            marker = list(color = c("#3498db", "#e74c3c", "#9b59b6"))) %>%
      layout(
        title = "Comparaison des MSE entre les estimateurs",
        yaxis = list(title = "MSE moyen"),
        xaxis = list(title = "Estimateur"),
        plot_bgcolor = "#f8f9fa"
      )
  })
  
  output$final_comparison_table <- DT::renderDataTable({
    df <- final_comparison_results()
    DT::datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 5))
  })
  
  
  optimal_params_stack <- reactiveValues(
    size = NULL,
    decay = NULL,
    maxit = NULL,
    h = NULL
  )
  
  optimal_params_stacking_nn <- reactiveValues(
    size = NA,
    decay = NA,
    maxit = NA
  )
  
  
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
  
  
  
  
  
  
  
  
  
  
  output$stack_plot_all <- renderPlotly({
    res <- stacking_results()
    if (is.null(res)) return()
    
    p <- plot_ly() %>%
      add_trace(x = res$x, y = res$m_true,
                type = 'scatter', mode = 'lines',
                name = "m(x) vraie",
                line = list(color = "black", dash = "dash", width = 2))
    
    method_colors <- c("nw" = "blue", "nn" = "red")
    legend_done <- c()
    
    for (method in res$methods_used) {
      if (!is.null(res$results[[method]])) {
        showlegend <- !(method %in% legend_done)  # Affiche la légende seulement une fois par méthode
        p <- p %>% add_trace(x = res$x, y = res$results[[method]],
                             type = 'scatter', mode = 'lines',
                             name = toupper(method),
                             line = list(color = method_colors[method]),
                             showlegend = showlegend)
        legend_done <- c(legend_done, method)
      }
    }
    
    # Ajout du stacking
    p <- p %>% add_trace(x = res$x, y = res$m_stack,
                         type = 'scatter', mode = 'lines',
                         name = paste("Stacking (", res$meta_method, ")"),
                         line = list(color = "green", width = 3))
    
    p %>% layout(
      title = "Comparaison des estimateurs avec Stacking",
      xaxis = list(title = "x"),
      yaxis = list(title = "m(x)"),
      legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center")
    )
  })
  
  
  output$stack_plot_facet <- renderPlotly({
    res <- stacking_results()
    if (is.null(res)) return()
    
    df_true <- data.frame(x = res$x, y = res$m_true, method = "m(x) vraie")
    p_list <- list()
    
    # Flag pour contrôler la légende
    legend_shown <- FALSE
    
    # NN vs m(x) vraie
    if ("nn" %in% res$methods_used && !is.null(res$results$nn)) {
      df_nn <- data.frame(x = res$x, y = res$results$nn, method = "NN")
      p_nn <- plot_ly() %>%
        add_trace(data = df_true, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
                  name = "m(x) vraie", showlegend = !legend_shown,
                  line = list(color = "black", dash = "dash")) %>%
        add_trace(data = df_nn, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
                  name = "NN", line = list(color = "red")) %>%
        layout(title = "NN vs m(x) vraie", yaxis = list(title = "m(x)"))
      p_list <- c(p_list, list(p_nn))
      legend_shown <- TRUE
    }
    
    # NW vs m(x) vraie
    if ("nw" %in% res$methods_used && !is.null(res$results$nw)) {
      df_nw <- data.frame(x = res$x, y = res$results$nw, method = "NW")
      p_nw <- plot_ly() %>%
        add_trace(data = df_true, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
                  name = "m(x) vraie", showlegend = !legend_shown,
                  line = list(color = "black", dash = "dash")) %>%
        add_trace(data = df_nw, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
                  name = "NW", line = list(color = "blue")) %>%
        layout(title = "NW vs m(x) vraie", yaxis = list(title = "m(x)"))
      p_list <- c(p_list, list(p_nw))
      legend_shown <- TRUE
    }
    
    # Stacking vs m(x) vraie
    if (!all(is.na(res$m_stack))) {
      stack_pred <- if (is.list(res$m_stack) || is.matrix(res$m_stack)) {
        as.numeric(res$m_stack)
      } else {
        res$m_stack
      }
      df_stack <- data.frame(x = res$x, y = stack_pred, method = paste("Stacking (", res$meta_method, ")"))
      p_stack <- plot_ly() %>%
        add_trace(data = df_true, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
                  name = "m(x) vraie", showlegend = !legend_shown,
                  line = list(color = "black", dash = "dash")) %>%
        add_trace(data = df_stack, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
                  name = "Stacking", line = list(color = "green")) %>%
        layout(title = paste("Stacking (", res$meta_method, ") vs m(x) vraie"), yaxis = list(title = "m(x)"))
      p_list <- c(p_list, list(p_stack))
      legend_shown <- TRUE
    }
    
    # Combinaison des facettes
    if (length(p_list) > 0) {
      subplot(p_list, nrows = length(p_list), shareX = TRUE, titleY = TRUE) %>%
        layout(showlegend = TRUE)
    }
  })
  
  
  optimal_params_stack <- reactiveValues(
    size = NULL,
    decay = NULL,
    maxit = NULL,
    h = NULL
  )
  
  
  stack_tuning_results <- eventReactive(input$run_stack_tuning, {
    req(input$stack_param, input$stack_nsim, input$stack_n, input$stacking_method, input$stack_methods)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    x_grid <- seq(0, 1, length.out = 100)
    h_grid <- seq(0.03, 0.4, length.out = 10)
    
    param_vals <- switch(input$stack_param,
                         "size" = c(1, 5, 10, 20, 30, 50),
                         "decay" = c(0.001, 0.01, 0.05, 0.1, 0.2),
                         "maxit" = c(100, 500, 1000, 3000))
    
    results <- list()
    h_used_vec <- numeric(input$stack_nsim)  # Pour mémoriser h utilisé si auto_h_stack
    
    withProgress(message = "Optimisation stacking...", value = 0, {
      for (i in seq_along(param_vals)) {
        pval <- param_vals[i]
        mse_vec <- numeric(input$stack_nsim)
        
        for (b in 1:input$stack_nsim) {
          data <- simulate_data(input$stack_n, m_true, input$t_dist, T_params)
          X_obs <- data$X_obs; Y_obs <- data$Y_obs; T_obs <- data$T_obs
          true_vals <- m_true(x_grid)
          
          df_meta <- data.frame(Y = Y_obs)
          m_nw_pred <- NULL; m_nn_pred <- NULL
          
          if ("nw" %in% input$stack_methods) {
            h_val <- if (isTRUE(input$auto_h_stack)) select_h_cv(X_obs, Y_obs, T_obs, h_grid) else input$stack_h
            h_used_vec[b] <- h_val  # Stocker h utilisé
            df_meta$m_nw <- nw_estimator_truncated(X_obs, X_obs, Y_obs, T_obs, h_val)
            m_nw_pred <- nw_estimator_truncated(x_grid, X_obs, Y_obs, T_obs, h_val)
          }
          
          if ("nn" %in% input$stack_methods) {
            nn_size <- if (input$stack_param == "size") pval else input$stack_nn_size
            nn_decay <- if (input$stack_param == "decay") pval else input$stack_nn_decay
            nn_maxit <- if (input$stack_param == "maxit") pval else input$stack_nn_maxit
            nn_res <- train_nn_estimator(X_obs, Y_obs, T_obs, size = nn_size, decay = nn_decay, maxit = nn_maxit, x_grid = x_grid)
            df_meta$m_nn <- predict(nn_res$model, newdata = data.frame(X = X_obs))
            m_nn_pred <- nn_res$predictions
          }
          
          pred_cols <- setdiff(names(df_meta), "Y")
          if (length(pred_cols) > 0) {
            df_pred <- data.frame(matrix(NA, nrow = length(x_grid), ncol = length(pred_cols)))
            colnames(df_pred) <- pred_cols
            for (col in pred_cols) {
              if (col == "m_nw") df_pred[[col]] <- m_nw_pred
              if (col == "m_nn") df_pred[[col]] <- m_nn_pred
            }
          } else {
            df_pred <- data.frame(dummy = rep(0, length(x_grid)))
          }
          
          meta_model <- switch(input$stacking_method,
                               "lm" = tryCatch(lm(Y ~ ., data = df_meta), error=function(e) NULL),
                               "mars" = tryCatch(earth(Y ~ ., data = df_meta), error=function(e) NULL),
                               "svr" = tryCatch(svm(Y ~ ., data = df_meta, kernel="radial", cost=1, gamma=1), error=function(e) NULL),
                               "xgb" = tryCatch({
                                 dtrain <- xgb.DMatrix(data = as.matrix(df_meta[, -1]), label = df_meta$Y)
                                 xgb.train(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose=0)
                               }, error=function(e) NULL),
                               "rf" = tryCatch(randomForest(Y ~ ., data = df_meta, ntree=500), error=function(e) NULL),
                               NULL)
          
          pred <- tryCatch({
            if (is.null(meta_model)) rep(NA, length(x_grid))
            else if (input$stacking_method == "xgb") predict(meta_model, newdata = as.matrix(df_pred))
            else predict(meta_model, newdata = df_pred)
          }, error = function(e) rep(NA, length(x_grid)))
          
          mse_vec[b] <- mean((pred - true_vals)^2, na.rm = TRUE)
        }
        
        results[[i]] <- data.frame(Param = pval, MSE = mean(mse_vec))
        incProgress(1/length(param_vals), detail = paste("Param =", pval))
      }
    })
    
    df_all <- do.call(rbind, results)
    # Filtrer les paramètres strictement positifs
    df_all <- df_all[df_all$Param > 0, ]
    
    if (nrow(df_all) == 0) {
      warning("Aucun paramètre strictement positif n'a été trouvé. Vérifie la définition de param_vals.")
      return(list(df = data.frame(), best_param = NA, best_mse = NA))
    }
    
    best_row <- df_all[which.min(df_all$MSE), ]
    
    if (input$stack_param == "size") {
      optimal_params_stack$size <- best_row$Param
    } else if (input$stack_param == "decay") {
      optimal_params_stack$decay <- best_row$Param
    } else if (input$stack_param == "maxit") {
      optimal_params_stack$maxit <- best_row$Param
    }
    if (isTRUE(input$auto_h_stack)) {
      optimal_params_stack$h <- ifelse(length(h_used_vec[h_used_vec > 0]) > 0, round(mean(h_used_vec[h_used_vec > 0], na.rm = TRUE), 3), NULL)
    }
    
    list(df = df_all, best_param = best_row$Param, best_mse = best_row$MSE)
  })
  
  
  
  
  
  
  
  
  
  
  
  optimal_params_stacking_nn <- reactiveValues(
    size = NA,
    decay = NA,
    maxit = NA
  )
  
  
  stacking_nn_param_results <- eventReactive(input$run_stacking_nn_study, {
    req(input$stacking_nn_param, input$stacking_nn_nsim)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    x_grid <- seq(0, 1, length.out = 100)
    
    param_vals <- switch(input$stacking_nn_param,
                         "size" = c(1, 5, 10, 20, 30, 50),
                         "decay" = c(0, 0.001, 0.01, 0.05, 0.1),
                         "maxit" = c(100, 500, 1000, 3000))
    
    grid_results <- list()
    
    withProgress(message = "Tuning NN pour Stacking...", value = 0, {
      for (i in seq_along(param_vals)) {
        pval <- param_vals[i]
        mse_vec <- numeric(input$stacking_nn_nsim)
        
        for (b in 1:input$stacking_nn_nsim) {
          data <- simulate_data(input$stacking_nn_n, m_true, input$t_dist, T_params)
          
          nn <- train_nn_estimator(
            X_obs = data$X_obs, Y_obs = data$Y_obs, T_obs = data$T_obs,
            size = if (input$stacking_nn_param == "size") pval else input$stacking_nn_fixed_size,
            decay = if (input$stacking_nn_param == "decay") pval else input$stacking_nn_fixed_decay,
            maxit = if (input$stacking_nn_param == "maxit") pval else input$stacking_nn_fixed_maxit,
            x_grid = x_grid
          )$predictions
          
          mse_vec[b] <- mean((nn - m_true(x_grid))^2, na.rm = TRUE)
        }
        
        grid_results[[i]] <- data.frame(Param = pval, MSE = mean(mse_vec))
        incProgress(1 / length(param_vals), detail = paste(input$stacking_nn_param, ":", pval))
      }
    })
    
    df_all <- do.call(rbind, grid_results)
    # Filtrer pour ne conserver que les paramètres strictement positifs
    df_all <- df_all[df_all$Param > 0, ]
    
    if (nrow(df_all) == 0) {
      warning("Aucun paramètre strictement positif trouvé. Vérifie la définition de param_vals.")
      return(list(df = data.frame(), best_param = NA, best_mse = NA))
    }
    
    best_row <- df_all[which.min(df_all$MSE), ]
    best_param <- best_row$Param
    
    optimal_params_stacking_nn[[input$stacking_nn_param]] <- best_param
    
    list(
      df = df_all,
      best_param = best_param,
      best_mse = best_row$MSE
    )
  })
  
  
  
  final_stack_results <- eventReactive(input$run_final_stack_comparison, {
    req(input$final_stack_n, input$final_stack_nsim, input$stack_methods, input$stacking_method)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    x_grid <- seq(0, 1, length.out = 100)
    h_grid <- seq(0.03, 0.4, length.out = 10)
    
    mse_nw <- numeric(input$final_stack_nsim)
    mse_nn <- numeric(input$final_stack_nsim)
    mse_stack <- numeric(input$final_stack_nsim)
    
    withProgress(message = "Simulation Stacking (NW vs NN vs Stacking)...", value = 0, {
      for (b in 1:input$final_stack_nsim) {
        data <- simulate_data(input$final_stack_n, m_true, input$t_dist, T_params)
        X_obs <- data$X_obs; Y_obs <- data$Y_obs; T_obs <- data$T_obs
        true_vals <- m_true(x_grid)
        
        # --- NW seul
        h_val <- if (isTRUE(input$final_auto_h)) {
          select_h_cv(X_obs, Y_obs, T_obs, h_grid)
        } else input$final_stack_h
        m_nw <- nw_estimator_truncated(x_grid, X_obs, Y_obs, T_obs, h_val)
        mse_nw[b] <- mean((m_nw - true_vals)^2, na.rm = TRUE)
        
        # --- NN seul
        nn_result <- train_nn_estimator(X_obs, Y_obs, T_obs,
                                        size = input$final_nn_size_solo, decay = input$final_nn_decay_solo,
                                        maxit = input$final_nn_maxit_solo, x_grid = x_grid)
        m_nn <- nn_result$predictions
        mse_nn[b] <- mean((m_nn - true_vals)^2, na.rm = TRUE)
        
        # --- Stacking (NN + NW)
        if (all(c("nw", "nn") %in% input$stack_methods)) {
          # Préparation meta-data
          df_meta <- data.frame(
            Y = Y_obs,
            m_nw = nw_estimator_truncated(X_obs, X_obs, Y_obs, T_obs, h_val),
            m_nn = predict(train_nn_estimator(X_obs, Y_obs, T_obs,
                                              size = input$final_nn_size_stack, decay = input$final_nn_decay_stack,
                                              maxit = input$final_nn_maxit_stack, x_grid = X_obs)$model,
                           newdata = data.frame(X = X_obs))
          )
          df_pred <- data.frame(
            m_nw = m_nw,
            m_nn = train_nn_estimator(X_obs, Y_obs, T_obs,
                                      size = input$final_nn_size_stack, decay = input$final_nn_decay_stack,
                                      maxit = input$final_nn_maxit_stack, x_grid = x_grid)$predictions
          )
          
          meta_model <- switch(input$stacking_method,
                               "lm" = tryCatch(lm(Y ~ ., data = df_meta), error=function(e) NULL),
                               "mars" = tryCatch(earth(Y ~ ., data = df_meta), error=function(e) NULL),
                               "svr" = tryCatch(svm(Y ~ ., data = df_meta, kernel="radial", cost=1, gamma=1), error=function(e) NULL),
                               "xgb" = tryCatch({
                                 dtrain <- xgb.DMatrix(data = as.matrix(df_meta[, -1]), label = df_meta$Y)
                                 xgb.train(data = dtrain, nrounds = input$final_xgb_nrounds, objective = "reg:squarederror", verbose=0)
                               }, error=function(e) NULL),
                               "rf" = tryCatch(randomForest(Y ~ ., data = df_meta, ntree=input$final_rf_ntree), error=function(e) NULL),
                               NULL)
          
          m_stack <- tryCatch({
            if (input$stacking_method == "xgb") {
              predict(meta_model, newdata = as.matrix(df_pred))
            } else {
              predict(meta_model, newdata = df_pred)
            }
          }, error=function(e) rep(NA, length(x_grid)))
          
          mse_stack[b] <- mean((m_stack - true_vals)^2, na.rm = TRUE)
        } else {
          mse_stack[b] <- NA
        }
        
        incProgress(1 / input$final_stack_nsim)
      }
    })
    
    data.frame(
      Méthode = c("NW seul", "NN seul", paste("Stacking (", input$stacking_method, ")", sep = "")),
      MSE_moyen = round(c(mean(mse_nw, na.rm=TRUE), mean(mse_nn, na.rm=TRUE), mean(mse_stack, na.rm=TRUE)), 5)
    )
  })
  
  output$final_stack_plot <- renderPlotly({
    res <- final_stack_results()
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    res <- res[!is.na(res$MSE_moyen), ]
    if (nrow(res) == 0) return(NULL)
    
    n_colors <- length(unique(res$Méthode))
    palette_colors <- if (n_colors == 1) c("blue") else if (n_colors == 2) c("blue", "orange") else RColorBrewer::brewer.pal(max(n_colors, 3), "Set2")
    
    plot_ly(data = res, x = ~Méthode, y = ~MSE_moyen, type = 'bar', color = ~Méthode,
            colors = palette_colors, text = ~paste("MSE moyen:", MSE_moyen), textposition = 'auto') %>%
      layout(title = "Comparaison MSE: NW seul vs NN seul vs Stacking",
             yaxis = list(rangemode = 'tozero', title = "MSE moyen"),
             xaxis = list(title = "Méthode"))
  })
  
  output$final_stack_table <- DT::renderDataTable({
    res <- final_stack_results()
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res <- res[!is.na(res$MSE_moyen), ]
    if (nrow(res) == 0) return(NULL)
    DT::datatable(res, rownames = FALSE, options = list(dom = 't', pageLength = 10)) %>%
      DT::formatRound('MSE_moyen', digits = 5)
  })
  
  
  
  
  
  
  
  
  
  
  
  
  # Résultats tuning du modèle méta (stacking NW + NN)
  # Résultats tuning du modèle méta (stacking)
  output$stack_param_plot <- renderPlotly({
    res <- stack_tuning_results()
    if (is.null(res)) return(NULL)
    
    plot_ly(res$df, x = ~Param, y = ~MSE, type = 'scatter', mode = 'lines+markers',
            name = "MSE", line = list(color = "purple")) %>%
      add_markers(x = res$best_param, y = res$best_mse, text = paste("Optimal:", res$best_param),
                  marker = list(color = "red", size = 10), name = "Optimal") %>%
      layout(title = paste("Optimisation du paramètre", input$stack_param),
             xaxis = list(title = input$stack_param),
             yaxis = list(title = "MSE moyen"),
             legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center"))
  })
  
  
  
  # Récapitulatif des paramètres optimaux du stacking
  output$summary_stacking_params <- DT::renderDataTable({
    df <- data.frame(
      Paramètre = c("size (neurones)", "decay (weight decay)", "maxit (itérations)", "h (bande de lissage)"),
      Valeur_optimale = c(
        ifelse(is.null(optimal_params_stack$size), "Non défini", round(optimal_params_stack$size, 3)),
        ifelse(is.null(optimal_params_stack$decay), "Non défini", round(optimal_params_stack$decay, 5)),
        ifelse(is.null(optimal_params_stack$maxit), "Non défini", round(optimal_params_stack$maxit, 0)),
        ifelse(is.null(optimal_params_stack$h), "Non défini", round(optimal_params_stack$h, 3))
      )
    )
    DT::datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 5))
  })
  
  # Meilleur paramètre sélectionné dans l’étude du stacking
  output$best_stack_param_display <- renderPrint({
    res <- stack_tuning_results()
    if (is.null(res)) {
      cat("Paramètre optimal non défini.")
    } else {
      cat("✓ Meilleur paramètre", input$stack_param, "=", round(res$best_param, 3),
          "\n✓ MSE moyen =", round(res$best_mse, 5))
    }
  })
  
  
  
  # Résumé textuel
  output$summary_stacking_text <- renderPrint({
    # Nettoyage et affichage unique des méthodes incluses
    methods_included <- if (!is.null(input$stack_methods)) unique(input$stack_methods) else "Aucune"
    methods_str <- paste(methods_included, collapse = ", ")
    
    cat("Résumé des paramètres optimaux du stacking :",
        "\n- Méthodes incluses :", methods_str,
        "\n- Méta-learner :", input$stacking_method)
  })
  
  
  # Paramètre optimal pour la composante NN dans le stacking
  output$best_stacking_nn_param_display <- renderPrint({
    res <- stacking_nn_param_results()
    if (is.null(res)) return(cat("Paramètre optimal NN non défini."))
    cat("✓ Meilleur paramètre", input$stacking_nn_param, "=", res$best_param)
  })
  
  # Plot pour tuning du réseau de neurones dans le stacking
  output$stacking_nn_param_plot <- renderPlotly({
    res <- stacking_nn_param_results()
    if (is.null(res)) return(NULL)
    
    df <- res$df
    
    plot_ly(df, x = ~Param) %>%
      add_trace(y = ~MSE, type = 'scatter', mode = 'lines+markers',
                name = "MSE", line = list(color = "purple")) %>%
      add_trace(
        x = res$best_param, y = res$best_mse, mode = "markers+text",
        text = paste0("Optimal = ", res$best_param),
        marker = list(color = "red", size = 10),
        name = "Optimal"
      ) %>%
      layout(
        title = paste("MSE NN selon", input$stacking_nn_param, "(Stacking)"),
        xaxis = list(title = input$stacking_nn_param),
        yaxis = list(title = "MSE moyen"),
        legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center")
      )
  })
  
  # Tableau des paramètres optimaux pour la partie NN dans le stacking
  output$summary_stacking_nn_params <- DT::renderDataTable({
    df <- data.frame(
      Paramètre = c("size (neurones)", "decay (weight decay)", "maxit (itérations)"),
      Valeur_optimale = c(
        ifelse(is.null(optimal_params_stacking_nn$size), "Non défini", optimal_params_stacking_nn$size),
        ifelse(is.null(optimal_params_stacking_nn$decay), "Non défini", optimal_params_stacking_nn$decay),
        ifelse(is.null(optimal_params_stacking_nn$maxit), "Non défini", optimal_params_stacking_nn$maxit)
      )
    )
    DT::datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 5))
  })
  
  # Confidence interval simulation
  ci_simulation <- eventReactive(input$run_ci_simulation, {
    req(input$ci_n, input$ci_b, input$ci_h, input$ci_level)
    
    m_true <- parse_true_func()
    T_params <- get_t_params()
    seed <- if (is.na(input$seed)) NULL else input$seed
    
    x_grid <- seq(0.05, 0.95, length.out = 100)
    z_alpha <- qnorm(1 - (1 - input$ci_level)/2)
    
    withProgress(message = 'Running CI simulation', value = 0, {
      m_hat_mat <- matrix(NA, nrow = input$ci_b, ncol = length(x_grid))
      
      for (b in 1:input$ci_b) {
        incProgress(1/input$ci_b, detail = paste("Simulation", b, "of", input$ci_b))
        if (!is.null(seed)) set.seed(seed + b)
        
        data <- simulate_data(input$ci_n, m_true, input$t_dist, T_params)
        if (input$auto_h_ci) {
          h_opt <- select_h_cv(data$X_obs, data$Y_obs, data$T_obs, h_grid)
        } else {
          h_opt <- input$ci_h
        }
        
        m_hat <- nw_estimator_truncated(x_grid, data$X_obs, data$Y_obs, data$T_obs, h_opt)
        m_hat_mat[b, ] <- m_hat
      }
    })
    
    # Calculate mean and CI
    m_hat_mean <- colMeans(m_hat_mat, na.rm = TRUE)
    m_hat_sd <- apply(m_hat_mat, 2, sd, na.rm = TRUE)
    ci_lower <- m_hat_mean - z_alpha * m_hat_sd
    ci_upper <- m_hat_mean + z_alpha * m_hat_sd
    
    # Calculate coverage
    true_vals <- m_true(x_grid)
    coverage <- mean(true_vals >= ci_lower & true_vals <= ci_upper, na.rm = TRUE)
    
    list(
      x = x_grid,
      m_true = true_vals,
      m_mean = m_hat_mean,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      coverage = coverage
    )
  })
  
  # Confidence interval plot
  output$ci_plot <- renderPlotly({
    results <- ci_simulation()
    if (is.null(results)) return()
    
    plot_ly() %>%
      add_trace(
        x = results$x, y = results$m_true,
        type = 'scatter', mode = 'lines',
        line = list(dash = "dash", color = "blue", width = 2),
        name = "True m(x)"
      ) %>%
      add_trace(
        x = results$x, y = results$m_mean,
        type = 'scatter', mode = 'lines',
        line = list(color = "black", width = 2),
        name = "Estimated mean"
      ) %>%
      add_ribbons(
        x = results$x,
        ymin = results$ci_lower,
        ymax = results$ci_upper,
        fillcolor = "rgba(255, 0, 0, 0.2)",
        line = list(color = "rgba(255, 0, 0, 0.1)"),
        name = paste0(input$ci_level*100, "% CI")
      ) %>%
      layout(
        title = paste("Confidence Intervals for NW Estimator with Truncation"),
        xaxis = list(title = "x"),
        yaxis = list(title = "m(x)"),
        legend = list(
          orientation = "h",
          x = 0.5,
          y = -0.3,
          xanchor = "center",
          yanchor = "top"
        ),
        margin = list(b = 100)
      )
  })
  
  output$ci_error_regions <- renderPlotly({
    results <- ci_simulation()
    if (is.null(results)) return()
    
    outside <- !(results$m_true >= results$ci_lower & results$m_true <= results$ci_upper)
    
    # Base plot
    p <- plot_ly() %>%
      add_trace(
        x = results$x, y = results$m_true,
        type = 'scatter', mode = 'lines',
        line = list(color = "blue", width = 2),
        name = "True m(x)"
      ) %>%
      add_ribbons(
        x = results$x,
        ymin = results$ci_lower,
        ymax = results$ci_upper,
        fillcolor = "rgba(255, 0, 0, 0.2)",
        line = list(color = "rgba(255, 0, 0, 0.1)"),
        name = paste0(input$ci_level * 100, "% CI")
      )
    
    # Ajouter uniquement s'il y a des points hors de l'IC
    if (any(outside, na.rm = TRUE)) {
      p <- p %>% add_markers(
        x = results$x[outside],
        y = results$m_true[outside],
        marker = list(color = "red", size = 8, symbol = "x"),
        name = "Outside CI"
      )
    }
    
    p %>% layout(
      title = "Zones de sous-couverture (Outside Confidence Intervals)",
      xaxis = list(title = "x"),
      yaxis = list(title = "m(x)"),
      legend = list(orientation = "h", x = 0.5, y = -0.3),
      margin = list(b = 100)
    )
  })
  
  
}

shinyApp(ui = ui, server = server)