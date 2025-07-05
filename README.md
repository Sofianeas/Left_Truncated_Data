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

To run the application directly from this repository, execute the following commands in RStudio:

```r
# In RStudio
if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
library(shiny)
runGitHub("Left_Truncated_Data", "Sofianeas", ref = "main")
```

### Running Locally

Alternatively, you can clone the repository and run it locally:

```bash
# in your terminal
git clone https://github.com/Sofianeas/Left_Truncated_Data.git
```

Then in R:

```r
# in R
if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
library(shiny)
setwd("path/to/Left_Truncated_Data")
runApp()
```

## 📊 Application Features

- **Estimation Panel**
  - Visualize the true function vs. estimated function  
  - Adjust bandwidth parameters  
  - Explore data before/after truncation  

- **Performance Analysis**
  - Study bias, variance, and MSE across sample sizes  
  - Compare estimator performance  

- **Bandwidth Study**
  - Examine how bandwidth selection affects estimation  
  - Find optimal bandwidth values  

- **MISE Simulation**
  - Analyze Mean Integrated Squared Error  
  - Study convergence properties  

- **Estimator Comparison**
  - Compare NW, LOESS, splines, and local polynomial estimators  
  - Visualize MSE differences  

- **NN + NW Combination**
  - Combine neural networks with kernel smoothing  
  - Optimize combination parameters  

- **Stacking Methods**
  - Implement meta-learning approaches  
  - Compare different stacking algorithms

## 📚 Theoretical Background

- Lynden–Bell estimator for truncation correction
- Nadaraya–Watson kernel regression adapted for truncated data
- Neural networks with truncation-aware weighting
- Stacking/ensemble methods for improved estimation



