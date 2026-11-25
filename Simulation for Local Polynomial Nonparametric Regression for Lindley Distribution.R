# SIMULATION STUDY
# Local Polynomial Nonparametric Regression under a
# Covariate-Dependent Lindley Distribution
#
# Local polynomial orders: p = 0, 1, 2
# Sample sizes: n = 50, 100, 200, 500
# Replications: R = 1000

rm(list = ls())

# 1. SIMULATION SETTINGS

set.seed(12345)

sample_sizes <- c(50, 100, 200, 500)

R <- 1000

degrees <- c(0, 1, 2)

# Constants in the true parameter function
a <- 0.5
b <- 0.5

# True parameter function
theta_true <- function(z) {
  exp(a + b * sin(2 * pi * z))
}

# Evaluation grid
G <- 51
z_grid <- seq(0, 1, length.out = G)

true_theta <- theta_true(z_grid)

# Candidate neighbourhood bandwidths
H <- seq(0.08, 0.50, length.out = 9)

# Numerical optimization settings
epsilon <- 1e-8
max_iterations <- 500


# 2. LINDLEY DISTRIBUTION

# Lindley probability density function
dLindley <- function(x, theta, log = FALSE) {
  
  log_density <-
    2 * log(theta) -
    log(theta + 1) +
    log1p(x) -
    theta * x
  
  if (log) {
    return(log_density)
  }
  
  return(exp(log_density))
}


# Generate random observations from the Lindley distribution
rLindley <- function(n, theta) {
  
  x <- numeric(n)
  
  # Mixture probability
  prob <- theta / (theta + 1)
  
  # Select exponential or gamma component
  component <- runif(n) < prob
  
  # Exponential component
  x[component] <- rexp(
    sum(component),
    rate = theta[component]
  )
  
  # Gamma component
  x[!component] <- rgamma(
    sum(!component),
    shape = 2,
    rate = theta[!component]
  )
  
  return(x)
}


# 3. GENERATE SIMULATED DATA

generate_data <- function(n) {
  
  # Generate covariates
  z <- runif(n, 0, 1)
  
  # Calculate true parameter values
  theta <- theta_true(z)
  
  # Generate Lindley responses
  x <- rLindley(n, theta)
  
  data.frame(
    X = x,
    Z = z,
    theta = theta
  )
}


# 4. LOCAL POLYNOMIAL DESIGN MATRIX

design_matrix <- function(z, z0, p) {
  
  d <- z - z0
  
  D <- sapply(0:p, function(j) d^j)
  
  D <- matrix(
    D,
    nrow = length(z),
    ncol = p + 1
  )
  
  return(D)
}


# 5. LOCAL POLYNOMIAL MAXIMUM LIKELIHOOD ESTIMATION

fit_local <- function(x, z, z0, h, p, exclude = NULL) {
  
  # Select observations inside the neighbourhood
  index <- which(abs(z - z0) <= h)
  
  # Exclude observation for LOOCV
  if (!is.null(exclude)) {
    index <- setdiff(index, exclude)
  }
  
  x_local <- x[index]
  z_local <- z[index]
  
  m <- length(index)
  
  # Check minimum local sample size
  if (m <= p + 1) {
    return(list(
      estimate = NA_real_,
      success = FALSE
    ))
  }
  
  # Construct local polynomial design matrix
  D <- design_matrix(z_local, z0, p)
  
  # Check rank of design matrix
  if (qr(D)$rank < p + 1) {
    return(list(
      estimate = NA_real_,
      success = FALSE
    ))
  }
  
  # Objective function: negative local log-likelihood
  objective <- function(beta) {
    
    theta <- as.vector(D %*% beta)
    
    # Positivity restriction
    if (any(theta <= 0) || any(!is.finite(theta))) {
      return(1e100)
    }
    
    loglik <- sum(
      2 * log(theta) -
        log(theta + 1) +
        log1p(x_local) -
        theta * x_local
    )
    
    return(-loglik)
  }
  
  # Analytical gradient of negative log-likelihood
  gradient <- function(beta) {
    
    theta <- as.vector(D %*% beta)
    
    if (any(theta <= 0) || any(!is.finite(theta))) {
      return(rep(0, p + 1))
    }
    
    g <- 2 / theta -
      1 / (theta + 1) -
      x_local
    
    grad <- -as.vector(crossprod(D, g))
    
    return(grad)
  }
  
  # Initial coefficient values
  beta_start <- c(
    max(1 / mean(x_local), 0.01),
    rep(0, p)
  )
  
  # Linear inequality constraints:
  # D %*% beta > epsilon
  ui <- D
  ci <- rep(epsilon, m)
  
  # Constrained numerical optimization
  fit <- tryCatch(
    constrOptim(
      theta = beta_start,
      f = objective,
      grad = gradient,
      ui = ui,
      ci = ci,
      method = "BFGS",
      control = list(
        maxit = max_iterations,
        reltol = 1e-9
      ),
      mu = 0.01
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(
      estimate = NA_real_,
      success = FALSE
    ))
  }
  
  # Check numerical convergence
  if (fit$convergence != 0 ||
      !is.finite(fit$value)) {
    
    return(list(
      estimate = NA_real_,
      success = FALSE
    ))
  }
  
  # Check positivity of fitted local parameters
  fitted_theta <- as.vector(D %*% fit$par)
  
  if (any(fitted_theta <= 0)) {
    return(list(
      estimate = NA_real_,
      success = FALSE
    ))
  }
  
  # Intercept estimates theta(z0)
  theta_hat <- fit$par[1]
  
  return(list(
    estimate = theta_hat,
    success = TRUE
  ))
}


# 6. LEAVE-ONE-OUT CROSS-VALIDATION

select_bandwidth <- function(x, z, p, H) {
  
  n <- length(x)
  
  cv_scores <- rep(Inf, length(H))
  
  for (k in seq_along(H)) {
    
    h <- H[k]
    
    cv_value <- 0
    valid <- TRUE
    
    for (i in seq_len(n)) {
      
      # Estimate theta(zi) after excluding observation i
      fit <- fit_local(
        x = x,
        z = z,
        z0 = z[i],
        h = h,
        p = p,
        exclude = i
      )
      
      if (!fit$success) {
        valid <- FALSE
        break
      }
      
      theta_hat <- fit$estimate
      
      # Negative log-likelihood of omitted observation
      contribution <- -dLindley(
        x[i],
        theta_hat,
        log = TRUE
      )
      
      if (!is.finite(contribution)) {
        valid <- FALSE
        break
      }
      
      cv_value <- cv_value + contribution
    }
    
    if (valid) {
      cv_scores[k] <- cv_value
    }
  }
  
  # Check whether any bandwidth gives valid CV fits
  if (all(!is.finite(cv_scores))) {
    
    return(list(
      h_opt = NA_real_,
      cv_scores = cv_scores,
      success = FALSE
    ))
  }
  
  # Select bandwidth with minimum CV score
  best <- which.min(cv_scores)
  
  return(list(
    h_opt = H[best],
    cv_scores = cv_scores,
    success = TRUE
  ))
}


# 7. STORAGE FOR SIMULATION RESULTS

simulation_results <- list()

summary_table <- data.frame()

bandwidth_table <- data.frame()

pointwise_table <- data.frame()

# 8. MAIN SIMULATION LOOP

for (n in sample_sizes) {
  
  cat("\n=\n")
  cat("Sample size:", n, "\n")
  cat("=\n")
  
  for (p in degrees) {
    
    cat("\nPolynomial degree:", p, "\n")
    
    # Store estimates for R replications and G grid points
    estimates <- matrix(
      NA_real_,
      nrow = R,
      ncol = G
    )
    
    # Store selected bandwidths
    selected_h <- rep(NA_real_, R)
    
    # Store CV success indicators
    cv_success <- logical(R)
    
   
    # REPLICATION LOOP
    
    for (r in seq_len(R)) {
      
      if (r %% 100 == 0) {
        cat("Replication:", r, "of", R, "\n")
      }
      
      # Generate simulated data
      dat <- generate_data(n)
      
      x <- dat$X
      z <- dat$Z
      
      # Select optimal bandwidth
      cv <- select_bandwidth(
        x = x,
        z = z,
        p = p,
        H = H
      )
      
      cv_success[r] <- cv$success
      
      if (!cv$success) {
        next
      }
      
      h_opt <- cv$h_opt
      
      selected_h[r] <- h_opt
      
      
      # ESTIMATE PARAMETER FUNCTION OVER THE GRID
      
      for (g in seq_len(G)) {
        
        fit <- fit_local(
          x = x,
          z = z,
          z0 = z_grid[g],
          h = h_opt,
          p = p
        )
        
        if (fit$success) {
          estimates[r, g] <- fit$estimate
        }
      }
    }
    
   
    # 9. CALCULATE PERFORMANCE MEASURES
    
    # Difference between estimated and true parameters
    error_matrix <- sweep(
      estimates,
      MARGIN = 2,
      STATS = true_theta,
      FUN = "-"
    )
    
    # Pointwise bias
    bias_grid <- colMeans(
      error_matrix,
      na.rm = TRUE
    )
    
    # Pointwise MSE
    mse_grid <- colMeans(
      error_matrix^2,
      na.rm = TRUE
    )
    
    # Pointwise MAE
    mae_grid <- colMeans(
      abs(error_matrix),
      na.rm = TRUE
    )
    
    # Replace undefined values if all fits failed at a grid point
    bias_grid[is.nan(bias_grid)] <- NA_real_
    mse_grid[is.nan(mse_grid)] <- NA_real_
    mae_grid[is.nan(mae_grid)] <- NA_real_
    
    # Integrated mean squared error
    imse <- mean(mse_grid, na.rm = TRUE)
    
    # Average pointwise bias
    average_bias <- mean(bias_grid, na.rm = TRUE)
    
    # Average absolute pointwise error
    average_mae <- mean(mae_grid, na.rm = TRUE)
    
    # Average selected bandwidth
    mean_h <- mean(selected_h, na.rm = TRUE)
    
    # Percentage of successful CV selections
    cv_rate <- mean(cv_success) * 100
    
    # Percentage of successful local fits
    convergence_rate <- mean(
      !is.na(estimates)
    ) * 100
    
    
    # 10. SAVE SUMMARY RESULTS
    
    summary_row <- data.frame(
      Sample_Size = n,
      Degree = p,
      Mean_Bandwidth = mean_h,
      Bias = average_bias,
      MAE = average_mae,
      IMSE = imse,
      CV_Success_Percent = cv_rate,
      Convergence_Percent = convergence_rate
    )
    
    summary_table <- rbind(
      summary_table,
      summary_row
    )
    
   
    # 11. SAVE POINTWISE PERFORMANCE
    
    pointwise_row <- data.frame(
      Sample_Size = n,
      Degree = p,
      Z = z_grid,
      True_Theta = true_theta,
      Bias = bias_grid,
      MSE = mse_grid,
      MAE = mae_grid
    )
    
    pointwise_table <- rbind(
      pointwise_table,
      pointwise_row
    )
    
    # 12. SAVE ESTIMATED CURVES
  
    simulation_results[[paste(n, p, sep = "_")]] <-
      list(
        estimates = estimates,
        mean_curve = colMeans(estimates, na.rm = TRUE),
        selected_h = selected_h,
        cv_success = cv_success,
        convergence_rate = convergence_rate
      )
    
    bandwidth_table <- rbind(
      bandwidth_table,
      data.frame(
        Sample_Size = n,
        Degree = p,
        Replication = seq_len(R),
        Selected_Bandwidth = selected_h,
        CV_Success = cv_success
      )
    )
  }
}


# 13. FINAL SIMULATION TABLE

summary_table <- summary_table[
  order(summary_table$Sample_Size,
        summary_table$Degree),
]

summary_table_rounded <- summary_table

numeric_columns <- vapply(
  summary_table_rounded,
  is.numeric,
  logical(1)
)

summary_table_rounded[numeric_columns] <-
  lapply(
    summary_table_rounded[numeric_columns],
    round,
    digits = 5
  )

print(summary_table_rounded)


# 14. EXPORT SIMULATION TABLES

write.csv(
  summary_table,
  "Simulation_Summary.csv",
  row.names = FALSE
)

write.csv(
  pointwise_table,
  "Pointwise_Performance.csv",
  row.names = FALSE
)

write.csv(
  bandwidth_table,
  "Selected_Bandwidths.csv",
  row.names = FALSE
)

# 15. GRAPHICAL REPRESENTATION
# True curve versus local polynomial estimates


# Four panels: n = 50, 100, 200, 500
# Each panel contains:
# True function, p = 0, p = 1, p = 2

png(
  filename = "Simulation_Parameter_Curves.png",
  width = 1400,
  height = 1100,
  res = 150
)

par(
  mfrow = c(2, 2),
  mar = c(4, 4, 3, 1)
)

for (n in sample_sizes) {
  
  key0 <- paste(n, 0, sep = "_")
  key1 <- paste(n, 1, sep = "_")
  key2 <- paste(n, 2, sep = "_")
  
  curve0 <- simulation_results[[key0]]$mean_curve
  curve1 <- simulation_results[[key1]]$mean_curve
  curve2 <- simulation_results[[key2]]$mean_curve
  
  # Determine plot limits
  all_values <- c(
    true_theta,
    curve0,
    curve1,
    curve2
  )
  
  y_limits <- range(
    all_values,
    finite = TRUE
  )
  
  plot(
    z_grid,
    true_theta,
    type = "l",
    lwd = 3,
    lty = 1,
    col = 1,
    ylim = y_limits,
    xlab = "Covariate (z)",
    ylab = expression(hat(theta)(z)),
    main = paste("Sample size n =", n)
  )
  
  lines(
    z_grid,
    curve0,
    col = 2,
    lwd = 2,
    lty = 2
  )
  
  lines(
    z_grid,
    curve1,
    col = 3,
    lwd = 2,
    lty = 3
  )
  
  lines(
    z_grid,
    curve2,
    col = 4,
    lwd = 2,
    lty = 4
  )
  
  legend(
    "topright",
    legend = c(
      "True parameter",
      "Local constant (p=0)",
      "Local linear (p=1)",
      "Local quadratic (p=2)"
    ),
    col = 1:4,
    lty = 1:4,
    lwd = c(3, 2, 2, 2),
    cex = 0.75,
    bty = "n"
  )
}

dev.off()

# 16. SAVE THE GRAPH IN PDF FORMAT

pdf(
  "Simulation_Parameter_Curves.pdf",
  width = 10,
  height = 8
)

par(
  mfrow = c(2, 2),
  mar = c(4, 4, 3, 1)
)

for (n in sample_sizes) {
  
  curve0 <- simulation_results[[
    paste(n, 0, sep = "_")
  ]]$mean_curve
  
  curve1 <- simulation_results[[
    paste(n, 1, sep = "_")
  ]]$mean_curve
  
  curve2 <- simulation_results[[
    paste(n, 2, sep = "_")
  ]]$mean_curve
  
  all_values <- c(
    true_theta,
    curve0,
    curve1,
    curve2
  )
  
  plot(
    z_grid,
    true_theta,
    type = "l",
    lwd = 3,
    col = 1,
    ylim = range(all_values, finite = TRUE),
    xlab = "Covariate (z)",
    ylab = expression(hat(theta)(z)),
    main = paste("Sample size n =", n)
  )
  
  lines(z_grid, curve0, col = 2, lwd = 2, lty = 2)
  lines(z_grid, curve1, col = 3, lwd = 2, lty = 3)
  lines(z_grid, curve2, col = 4, lwd = 2, lty = 4)
  
  legend(
    "topright",
    legend = c(
      "True parameter",
      "Local constant (p=0)",
      "Local linear (p=1)",
      "Local quadratic (p=2)"
    ),
    col = 1:4,
    lty = 1:4,
    lwd = c(3, 2, 2, 2),
    cex = 0.75,
    bty = "n"
  )
}

dev.off()

# 17. PLOT THE TRUE PARAMETER FUNCTION ALONE

png(
  "True_Parameter_Function.png",
  width = 1000,
  height = 700,
  res = 150
)

plot(
  z_grid,
  true_theta,
  type = "l",
  lwd = 3,
  xlab = "Covariate (z)",
  ylab = expression(theta(z)),
  main = expression(
    theta(z) == exp(0.5 + 0.5 * sin(2 * pi * z))
  )
)

dev.off()

# 18. PRINT FINAL RESULTS

cat("\n=\n")
cat("SIMULATION STUDY COMPLETED\n")
cat("=\n")

cat("\nTrue parameter function:\n")
cat("theta(z) = exp(0.5 + 0.5*sin(2*pi*z))\n")

cat("\nNumber of replications:", R, "\n")

cat("\nSimulation summary:\n")
print(summary_table_rounded)

cat("\nFiles generated:\n")
cat("1. Simulation_Summary.csv\n")
cat("2. Pointwise_Performance.csv\n")
cat("3. Selected_Bandwidths.csv\n")
cat("4. Simulation_Parameter_Curves.png\n")
cat("5. Simulation_Parameter_Curves.pdf\n")
cat("6. True_Parameter_Function.png\n")


# 19. ADDITIONAL GRAPHS FOR SIMULATION RESULTS

# Install ggplot2 and tidyr if required
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

library(ggplot2)
library(tidyr)

# Read actual simulation results
results <- read.csv("Simulation_Summary.csv")

pointwise <- read.csv("Pointwise_Performance.csv")

bandwidth_data <- read.csv("Selected_Bandwidths.csv")

# Convert polynomial degree into a factor
results$Degree <- factor(
  results$Degree,
  levels = c(0, 1, 2),
  labels = c("p = 0", "p = 1", "p = 2")
)

pointwise$Degree <- factor(
  pointwise$Degree,
  levels = c(0, 1, 2),
  labels = c("p = 0", "p = 1", "p = 2")
)

bandwidth_data$Degree <- factor(
  bandwidth_data$Degree,
  levels = c(0, 1, 2),
  labels = c("p = 0", "p = 1", "p = 2")
)

# Common theme
my_theme <- theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 13
    ),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# GRAPH 1: BIAS COMPARISON

g1 <- ggplot(
  results,
  aes(
    x = Sample_Size,
    y = Bias,
    color = Degree,
    group = Degree
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = sort(unique(results$Sample_Size))
  ) +
  labs(
    title = "Comparison of Bias",
    x = "Sample Size (n)",
    y = "Average Bias"
  ) +
  my_theme

print(g1)

ggsave(
  "Graph1_Bias.png",
  g1,
  width = 7,
  height = 5,
  dpi = 300
)


# GRAPH 2: MAE COMPARISON

g2 <- ggplot(
  results,
  aes(
    x = Sample_Size,
    y = MAE,
    color = Degree,
    group = Degree
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = sort(unique(results$Sample_Size))
  ) +
  labs(
    title = "Comparison of Mean Absolute Error",
    x = "Sample Size (n)",
    y = "MAE"
  ) +
  my_theme

print(g2)

ggsave(
  "Graph2_MAE.png",
  g2,
  width = 7,
  height = 5,
  dpi = 300
)

# GRAPH 3: IMSE COMPARISON

g3 <- ggplot(
  results,
  aes(
    x = Sample_Size,
    y = IMSE,
    color = Degree,
    group = Degree
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = sort(unique(results$Sample_Size))
  ) +
  labs(
    title = "Comparison of Integrated Mean Squared Error",
    x = "Sample Size (n)",
    y = "IMSE"
  ) +
  my_theme

print(g3)

ggsave(
  "Graph3_IMSE.png",
  g3,
  width = 7,
  height = 5,
  dpi = 300
)

# GRAPH 4: MEAN BANDWIDTH COMPARISON

g4 <- ggplot(
  results,
  aes(
    x = Sample_Size,
    y = Mean_Bandwidth,
    color = Degree,
    group = Degree
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = sort(unique(results$Sample_Size))
  ) +
  labs(
    title = "Mean Selected Bandwidth",
    x = "Sample Size (n)",
    y = "Mean Bandwidth"
  ) +
  my_theme

print(g4)

ggsave(
  "Graph4_Bandwidth.png",
  g4,
  width = 7,
  height = 5,
  dpi = 300
)

# GRAPH 5: CV SUCCESS AND CONVERGENCE PERCENTAGES

reliability_data <- results |>
  pivot_longer(
    cols = c(
      CV_Success_Percent,
      Convergence_Percent
    ),
    names_to = "Measure",
    values_to = "Percentage"
  )

reliability_data$Measure <- factor(
  reliability_data$Measure,
  levels = c(
    "CV_Success_Percent",
    "Convergence_Percent"
  ),
  labels = c(
    "CV Success",
    "Convergence"
  )
)

g5 <- ggplot(
  reliability_data,
  aes(
    x = Sample_Size,
    y = Percentage,
    color = Degree,
    group = Degree
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~Measure) +
  scale_x_continuous(
    breaks = sort(unique(results$Sample_Size))
  ) +
  labs(
    title = "Cross-Validation Success and Convergence",
    x = "Sample Size (n)",
    y = "Percentage (%)"
  ) +
  my_theme

print(g5)

ggsave(
  "Graph5_CV_Convergence.png",
  g5,
  width = 9,
  height = 5,
  dpi = 300
)


# GRAPH 6: POINTWISE BIAS CURVES

g6 <- ggplot(
  pointwise,
  aes(
    x = Z,
    y = Bias,
    color = Degree
  )
) +
  geom_line(linewidth = 0.9) +
  facet_wrap(
    ~Sample_Size,
    scales = "free_y"
  ) +
  labs(
    title = "Pointwise Bias Curves",
    x = "Covariate (z)",
    y = "Pointwise Bias"
  ) +
  my_theme

print(g6)

ggsave(
  "Graph6_Pointwise_Bias.png",
  g6,
  width = 10,
  height = 7,
  dpi = 300
)


# GRAPH 7: POINTWISE MSE CURVES

g7 <- ggplot(
  pointwise,
  aes(
    x = Z,
    y = MSE,
    color = Degree
  )
) +
  geom_line(linewidth = 0.9) +
  facet_wrap(
    ~Sample_Size,
    scales = "free_y"
  ) +
  labs(
    title = "Pointwise Mean Squared Error",
    x = "Covariate (z)",
    y = "Pointwise MSE"
  ) +
  my_theme

print(g7)

ggsave(
  "Graph7_Pointwise_MSE.png",
  g7,
  width = 10,
  height = 7,
  dpi = 300
)

# GRAPH 8: POINTWISE MAE CURVES

g8 <- ggplot(
  pointwise,
  aes(
    x = Z,
    y = MAE,
    color = Degree
  )
) +
  geom_line(linewidth = 0.9) +
  facet_wrap(
    ~Sample_Size,
    scales = "free_y"
  ) +
  labs(
    title = "Pointwise Mean Absolute Error",
    x = "Covariate (z)",
    y = "Pointwise MAE"
  ) +
  my_theme

print(g8)

ggsave(
  "Graph8_Pointwise_MAE.png",
  g8,
  width = 10,
  height = 7,
  dpi = 300
)

# GRAPH 9: DISTRIBUTION OF SELECTED BANDWIDTHS

g9 <- ggplot(
  bandwidth_data,
  aes(
    x = Degree,
    y = Selected_Bandwidth,
    fill = Degree
  )
) +
  geom_boxplot(
    na.rm = TRUE,
    alpha = 0.7
  ) +
  facet_wrap(~Sample_Size) +
  labs(
    title = "Distribution of Selected Bandwidths",
    x = "Polynomial Degree",
    y = "Selected Bandwidth"
  ) +
  my_theme +
  theme(legend.position = "none")

print(g9)

ggsave(
  "Graph9_Bandwidth_Distribution.png",
  g9,
  width = 10,
  height = 7,
  dpi = 300
)

# END OF ADDITIONAL GRAPHS

cat("\nAll additional simulation graphs have been generated.\n")

