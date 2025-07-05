# Nonparametric Estimation with Left-Truncated Data – Shiny Application
This application allows you to simulate left-truncated data, estimate the associated regression function and compare between several estimators
## 📌 Overview
This repository contains an interactive Shiny application for nonparametric estimation with left-truncated data. The app provides a comprehensive toolkit for:

- Estimating regression functions with left-truncated observations  
- Comparing various nonparametric estimators  
- Visualizing the impact of bandwidth selection  
- Analyzing estimator performance through simulations  
- Combining estimators (Nadaraya–Watson and neural networks) via stacking  

## 🚀 How to Run the App

> **Important Notice Before Running**  
> This application requires significant computational resources for some simulations. For optimal performance:
> - Start with smaller sample sizes (n < 500) when first testing the app  
> - Reduce the number of simulations (N_sim < 50) for performance analysis  
> - Be patient – some calculations (especially MISE and uniform convergence) may take several minutes  
> - Monitor your system resources – the app may use substantial memory for large simulations  

### Running from GitHub

```r
# In RStudio
if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
library(shiny)
runGitHub("Left_Truncated_Data", "Sofianeas", ref = "main")

