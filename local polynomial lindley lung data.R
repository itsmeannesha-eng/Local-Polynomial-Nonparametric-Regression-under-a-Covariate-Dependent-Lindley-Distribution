#Complete R Script: Real Data Analysis & Model Comparison on the `lung` Dataset
# Paper: Local Polynomial Nonparametric Regression under a Covariate-Dependent
#        Lindley Distribution

# 1. Load Required Packages & Clean Data (Sections 2.10 & 5.1)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)

data(lung, package = "survival")

df <- lung %>%
  select(time, age) %>%
  filter(!is.na(time), !is.na(age), is.finite(time), is.finite(age), time > 0)

# Scale survival time to months for numerical stability of rate parameter theta
df$X <- df$time / 30.4375   # Strictly positive response: Survival time (months)
df$Z <- df$age              # Continuous covariate: Age (years)
n <- nrow(df)

# 2. Core Functions: Lindley Density, Local Polynomial MLE & LOO-CV (Section 2)

# Lindley log-density: Eq. (2.15)
lindley_logpdf <- function(x, theta) {
  2 * log(theta) - log(theta + 1) + log(1 + x) - theta * x
}

# Conditional mean E(X | Z = z): Eq. (2.6)
lindley_mean <- function(theta) {
  (theta + 2) / (theta * (theta + 1))
}

# Local polynomial MLE at target point z0 with bandwidth h: Eq. (2.16)-(2.18)
fit_local_lindley_point <- function(x, z, z0, h, p = 1) {
  idx <- which(abs(z - z0) <= h)
  n_loc <- length(idx)
  
  # Require local sample size > p + 1 (Section 2.7)
  if (n_loc <= (p + 1)) {
    return(list(theta0 = NA, beta = rep(NA, p + 1), converged = FALSE, n_loc = n_loc, loglik = NA))
  }
  
  d <- z[idx] - z0
  x_loc <- x[idx]
  
  D_mat <- sapply(0:p, function(j) d^j)
  if (!is.matrix(D_mat)) D_mat <- matrix(D_mat, ncol = 1)
  
  if (qr(D_mat)$rank < (p + 1)) {
    return(list(theta0 = NA, beta = rep(NA, p + 1), converged = FALSE, n_loc = n_loc, loglik = NA))
  }
  
  # Moment-based initial guess for intercept
  m_loc <- mean(x_loc)
  theta_init <- max(1e-3, (-(m_loc - 1) + sqrt((m_loc - 1)^2 + 8 * m_loc)) / (2 * m_loc))
  beta_init <- c(theta_init, rep(0, p))
  
  neg_ll <- function(beta) {
    theta_vec <- as.vector(D_mat %*% beta)
    if (any(!is.finite(theta_vec)) || any(theta_vec <= 1e-6)) {
      return(1e12)
    }
    -sum(lindley_logpdf(x_loc, theta_vec))
  }
  
  opt <- optim(
    par = beta_init,
    fn = neg_ll,
    method = "BFGS",
    control = list(maxit = 500, reltol = 1e-8)
  )
  
  theta_hat <- as.vector(D_mat %*% opt$par)
  valid <- (opt$convergence == 0) && all(theta_hat > 1e-6) && is.finite(opt$value)
  
  list(
    theta0 = ifelse(valid, opt$par[1], NA),
    beta = if (valid) opt$par else rep(NA, p + 1),
    converged = valid,
    n_loc = n_loc,
    loglik = ifelse(valid, -opt$value, NA)
  )
}

# Leave-One-Out Cross-Validation: Eq. (2.26)-(2.27)
cv_bandwidth_select <- function(x, z, h_grid, p = 1) {
  n_obs <- length(x)
  res <- data.frame(h = h_grid, CV = NA, valid_rate = NA)
  
  for (k in seq_along(h_grid)) {
    h <- h_grid[k]
    nll_loo <- numeric(n_obs)
    ok <- logical(n_obs)
    
    for (i in 1:n_obs) {
      fit_i <- fit_local_lindley_point(x[-i], z[-i], z0 = z[i], h = h, p = p)
      if (fit_i$converged && !is.na(fit_i$theta0)) {
        nll_loo[i] <- -lindley_logpdf(x[i], fit_i$theta0)
        ok[i] <- TRUE
      } else {
        ok[i] <- FALSE
      }
    }
    
    res$valid_rate[k] <- mean(ok) * 100
    if (all(ok)) {
      res$CV[k] <- sum(nll_loo)
    }
  }
  res
}

# Evaluate fitted curve across a grid of target covariate values: Eq. (2.25)
fit_curve <- function(x, z, z_grid, h, p = 1) {
  fits <- lapply(z_grid, function(z0) fit_local_lindley_point(x, z, z0, h, p))
  data.frame(
    z = z_grid,
    theta_hat = sapply(fits, `[[`, "theta0"),
    mean_hat = lindley_mean(sapply(fits, `[[`, "theta0")),
    converged = sapply(fits, `[[`, "converged"),
    n_loc = sapply(fits, `[[`, "n_loc"),
    degree = paste0("p = ", p)
  )
}

# 3. Table 3: Descriptive Statistics (Section 5.1)

desc_stats <- function(v, name) {
  data.frame(
    Variable = name,
    n = length(v),
    Minimum = round(min(v), 4),
    Maximum = round(max(v), 4),
    Mean = round(mean(v), 4),
    Median = round(median(v), 4),
    SD = round(sd(v), 4)
  )
}

table3_descriptive <- rbind(
  desc_stats(df$X, "Response X (Survival Time, Months)"),
  desc_stats(df$Z, "Covariate Z (Age, Years)")
)

cat("\n--- TABLE 3: DESCRIPTIVE STATISTICS ---\n")
print(table3_descriptive)

# 4. Bandwidth Selection & Model Comparison Table (Sections 2.9 & 5.2)

h_grid <- seq(6, 22, by = 2)

cv_p0 <- cv_bandwidth_select(df$X, df$Z, h_grid, p = 0)
cv_p1 <- cv_bandwidth_select(df$X, df$Z, h_grid, p = 1)
cv_p2 <- cv_bandwidth_select(df$X, df$Z, h_grid, p = 2)

h_opt_0 <- cv_p0$h[which.min(cv_p0$CV)]
h_opt_1 <- cv_p1$h[which.min(cv_p1$CV)]
h_opt_2 <- cv_p2$h[which.min(cv_p2$CV)]

eval_local_model <- function(x, z, h, p, label) {
  n_obs <- length(x)
  theta_in <- numeric(n_obs)
  nll_loo <- numeric(n_obs)
  conv <- logical(n_obs)
  
  for (i in 1:n_obs) {
    f_in <- fit_local_lindley_point(x, z, z[i], h, p)
    f_loo <- fit_local_lindley_point(x[-i], z[-i], z[i], h, p)
    theta_in[i] <- f_in$theta0
    conv[i] <- f_in$converged
    nll_loo[i] <- -lindley_logpdf(x[i], f_loo$theta0)
  }
  
  ll <- sum(lindley_logpdf(x, theta_in))
  edf <- (p + 1) * (diff(range(z)) / (2 * h))
  aic <- -2 * ll + 2 * edf
  bic <- -2 * ll + log(n_obs) * edf
  pred_mean <- lindley_mean(theta_in)
  
  data.frame(
    Model = label,
    Bandwidth_or_Params = paste0("h = ", h, " (edf=", round(edf, 2), ")"),
    LogLik = round(ll, 3),
    AIC = round(aic, 3),
    BIC = round(bic, 3),
    LOO_NLL = round(sum(nll_loo), 3),
    RMSE = round(sqrt(mean((x - pred_mean)^2)), 4),
    MAE = round(mean(abs(x - pred_mean)), 4),
    Convergence_Pct = round(mean(conv) * 100, 2)
  )
}

# Competing Parametric Models
# (a) Constant Lindley
opt_const <- optimize(function(th) -sum(lindley_logpdf(df$X, th)), interval = c(1e-4, 10))
ll_const <- -opt_const$objective
mean_const <- lindley_mean(opt_const$minimum)
loo_const <- sum(sapply(1:n, function(i) {
  th_i <- optimize(function(th) -sum(lindley_logpdf(df$X[-i], th)), interval = c(1e-4, 10))$minimum
  -lindley_logpdf(df$X[i], th_i)
}))

# (b) Parametric Log-Linear Lindley: theta(z) = exp(b0 + b1*z)
opt_param_lind <- optim(c(log(opt_const$minimum), 0), function(b) {
  -sum(lindley_logpdf(df$X, exp(b[1] + b[2] * df$Z)))
})
ll_pl <- -opt_param_lind$value
th_pl <- exp(opt_param_lind$par[1] + opt_param_lind$par[2] * df$Z)
mean_pl <- lindley_mean(th_pl)
loo_pl <- sum(sapply(1:n, function(i) {
  b_i <- optim(opt_param_lind$par, function(b) -sum(lindley_logpdf(df$X[-i], exp(b[1] + b[2] * df$Z[-i]))))$par
  -lindley_logpdf(df$X[i], exp(b_i[1] + b_i[2] * df$Z[i]))
}))

# (c) Exponential Regression
fit_exp <- survreg(Surv(X) ~ Z, data = df, dist = "exponential")
ll_exp <- as.numeric(logLik(fit_exp))
mean_exp <- predict(fit_exp, type = "response")
loo_exp <- sum(sapply(1:n, function(i) {
  m_i <- survreg(Surv(X) ~ Z, data = df[-i, ], dist = "exponential")
  mu_i <- predict(m_i, newdata = df[i, ], type = "response")
  -dexp(df$X[i], rate = 1 / mu_i, log = TRUE)
}))

# (d) Gamma Regression (log link)
fit_gam <- glm(X ~ Z, data = df, family = Gamma(link = "log"))
ll_gam <- as.numeric(logLik(fit_gam))
mean_gam <- predict(fit_gam, type = "response")
loo_gam <- sum(sapply(1:n, function(i) {
  m_i <- glm(X ~ Z, data = df[-i, ], family = Gamma(link = "log"))
  mu_i <- predict(m_i, newdata = df[i, ], type = "response")
  sh_i <- 1 / summary(m_i)$dispersion
  -dgamma(df$X[i], shape = sh_i, scale = mu_i / sh_i, log = TRUE)
}))

# (e) Weibull Regression
fit_wei <- survreg(Surv(X) ~ Z, data = df, dist = "weibull")
ll_wei <- as.numeric(logLik(fit_wei))
mean_wei <- predict(fit_wei, type = "response") * gamma(1 + fit_wei$scale)
loo_wei <- sum(sapply(1:n, function(i) {
  m_i <- survreg(Surv(X) ~ Z, data = df[-i, ], dist = "weibull")
  sc_i <- predict(m_i, newdata = df[i, ], type = "response")
  sh_i <- 1 / m_i$scale
  -dweibull(df$X[i], shape = sh_i, scale = sc_i, log = TRUE)
}))

table_model_comparison <- rbind(
  eval_local_model(df$X, df$Z, h_opt_0, p = 0, "Local Constant Lindley (p = 0)"),
  eval_local_model(df$X, df$Z, h_opt_1, p = 1, "Local Linear Lindley (p = 1)"),
  eval_local_model(df$X, df$Z, h_opt_2, p = 2, "Local Quadratic Lindley (p = 2)"),
  data.frame(Model = "Constant Lindley (Parametric)", Bandwidth_or_Params = "1 param",
             LogLik = round(ll_const, 3), AIC = round(-2*ll_const + 2, 3), BIC = round(-2*ll_const + log(n), 3),
             LOO_NLL = round(loo_const, 3), RMSE = round(sqrt(mean((df$X - mean_const)^2)), 4),
             MAE = round(mean(abs(df$X - mean_const)), 4), Convergence_Pct = 100),
  data.frame(Model = "Log-Linear Lindley (Parametric)", Bandwidth_or_Params = "2 params",
             LogLik = round(ll_pl, 3), AIC = round(-2*ll_pl + 4, 3), BIC = round(-2*ll_pl + 2*log(n), 3),
             LOO_NLL = round(loo_pl, 3), RMSE = round(sqrt(mean((df$X - mean_pl)^2)), 4),
             MAE = round(mean(abs(df$X - mean_pl)), 4), Convergence_Pct = 100),
  data.frame(Model = "Exponential Regression", Bandwidth_or_Params = "2 params",
             LogLik = round(ll_exp, 3), AIC = round(AIC(fit_exp), 3), BIC = round(BIC(fit_exp), 3),
             LOO_NLL = round(loo_exp, 3), RMSE = round(sqrt(mean((df$X - mean_exp)^2)), 4),
             MAE = round(mean(abs(df$X - mean_exp)), 4), Convergence_Pct = 100),
  data.frame(Model = "Gamma Regression", Bandwidth_or_Params = "3 params",
             LogLik = round(ll_gam, 3), AIC = round(AIC(fit_gam), 3), BIC = round(BIC(fit_gam), 3),
             LOO_NLL = round(loo_gam, 3), RMSE = round(sqrt(mean((df$X - mean_gam)^2)), 4),
             MAE = round(mean(abs(df$X - mean_gam)), 4), Convergence_Pct = 100),
  data.frame(Model = "Weibull Regression", Bandwidth_or_Params = "3 params",
             LogLik = round(ll_wei, 3), AIC = round(AIC(fit_wei), 3), BIC = round(BIC(fit_wei), 3),
             LOO_NLL = round(loo_wei, 3), RMSE = round(sqrt(mean((df$X - mean_wei)^2)), 4),
             MAE = round(mean(abs(df$X - mean_wei)), 4), Convergence_Pct = 100)
)

cat("\n--- MODEL COMPARISON TABLE ---\n")
print(table_model_comparison)

# 5. Table 2: Sensitivity Analysis Summary (Section 4)

z_grid <- seq(min(df$Z), max(df$Z), length.out = 100)

sens_bw <- rbind(
  cbind(Setting = "Small Bandwidth (h = 6)", fit_curve(df$X, df$Z, z_grid, h = 6, p = 1)),
  cbind(Setting = paste0("Optimal Bandwidth (h = ", h_opt_1, ")"), fit_curve(df$X, df$Z, z_grid, h = h_opt_1, p = 1)),
  cbind(Setting = "Large Bandwidth (h = 20)", fit_curve(df$X, df$Z, z_grid, h = 20, p = 1))
)

table2_sensitivity <- sens_bw %>%
  group_by(Setting) %>%
  summarise(
    Min_Theta = round(min(theta_hat, na.rm = TRUE), 4),
    Max_Theta = round(max(theta_hat, na.rm = TRUE), 4),
    Curve_SD = round(sd(theta_hat, na.rm = TRUE), 4),
    Min_Local_n = min(n_loc),
    Median_Local_n = median(n_loc),
    Convergence_Pct = round(mean(converged) * 100, 2)
  )

cat("\n--- TABLE 2: SENSITIVITY ANALYSIS ---\n")
print(table2_sensitivity)

# 6. All Diagnostic & Comparison Figures (Sections 4 & 5)

# Plot 1: Preliminary Scatter Plot (Section 5.1)
p1_scatter <- ggplot(df, aes(x = Z, y = X)) +
  geom_point(alpha = 0.65, color = "#2C3E50", size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "#E74C3C", fill = "#FADBD8") +
  labs(
    title = "",
    x = "Age (Years)",
    y = "Survival Time (Months)"
  ) +
  theme_bw(base_size = 12)

# Plot 2 (Figure 2 in paper): Estimated Parameter Function theta_hat(z) across p = 0, 1, 2
curve_p0 <- fit_curve(df$X, df$Z, z_grid, h = h_opt_0, p = 0)
curve_p1 <- fit_curve(df$X, df$Z, z_grid, h = h_opt_1, p = 1)
curve_p2 <- fit_curve(df$X, df$Z, z_grid, h = h_opt_2, p = 2)
all_curves <- rbind(curve_p0, curve_p1, curve_p2)

p2_theta_curve <- ggplot(all_curves, aes(x = z, y = theta_hat, color = degree, linetype = degree)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = expression(paste("", )),
    x = "Age (Years)",
    y = expression(hat(theta)(z)),
    color = "Polynomial Degree",
    linetype = "Polynomial Degree"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Plot 3: Conditional Mean Survival E(X | Z = z) Comparison Across Models
cond_mean_df <- rbind(
  data.frame(z = curve_p1$z, mean_val = curve_p1$mean_hat, Model = "Local Linear Lindley (p = 1)"),
  data.frame(z = curve_p2$z, mean_val = curve_p2$mean_hat, Model = "Local Quadratic Lindley (p = 2)"),
  data.frame(z = z_grid, mean_val = lindley_mean(exp(opt_param_lind$par[1] + opt_param_lind$par[2] * z_grid)),
             Model = "Log-Linear Lindley (Parametric)"),
  data.frame(z = z_grid, mean_val = predict(fit_gam, newdata = data.frame(Z = z_grid), type = "response"),
             Model = "Gamma Regression"),
  data.frame(z = z_grid, mean_val = predict(fit_wei, newdata = data.frame(Z = z_grid), type = "response") * gamma(1 + fit_wei$scale),
             Model = "Weibull Regression")
)

p3_cond_mean <- ggplot() +
  geom_point(data = df, aes(x = Z, y = X), alpha = 0.25, color = "gray40") +
  geom_line(data = cond_mean_df, aes(x = z, y = mean_val, color = Model), linewidth = 1.05) +
  labs(
    title = "",
    x = "Age (Years)",
    y = "Estimated Mean Survival Time (Months)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Plot 4: Leave-One-Out Cross-Validation Criterion vs. Bandwidth (Section 2.9)
cv_all <- rbind(
  cbind(cv_p0, Degree = "p = 0"),
  cbind(cv_p1, Degree = "p = 1"),
  cbind(cv_p2, Degree = "p = 2")
) %>% filter(!is.na(CV))

p4_cv_curve <- ggplot(cv_all, aes(x = h, y = CV, color = Degree)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  labs(
    title = "",
    x = "Bandwidth (h)",
    y = "LOO Negative Log-Likelihood"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Plot 5: Sensitivity of Estimated Parameter Curve to Bandwidth Selection (Section 4.2)
p5_sens_bw <- ggplot(sens_bw, aes(x = z, y = theta_hat, color = Setting)) +
  geom_line(linewidth = 1.05) +
  labs(
    title = "",
    x = "Age (Years)",
    y = expression(hat(theta)(z))
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Plot 6: Boundary Behaviour & Local Sample Size Diagnostic (Section 4.4)
p6_boundary <- ggplot(curve_p1, aes(x = z, y = n_loc)) +
  geom_line(color = "#2980B9", linewidth = 1.1) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  labs(
    title = paste0(""),
    x = "Age (Years)",
    y = "Number of Observations in N(z0, h)"
  ) +
  theme_bw(base_size = 12)

# Plot 7: Conditional Lindley Densities at Representative Ages (Section 5.3)
ages_rep <- quantile(df$Z, probs = c(0.15, 0.50, 0.85))
x_vals <- seq(0.01, max(df$X), length.out = 200)

dens_df <- do.call(rbind, lapply(ages_rep, function(a) {
  th <- fit_local_lindley_point(df$X, df$Z, z0 = a, h = h_opt_1, p = 1)$theta0
  data.frame(
    x = x_vals,
    density = exp(lindley_logpdf(x_vals, th)),
    Age_Group = paste0("Age = ", round(a, 1), " yrs (theta = ", round(th, 3), ")")
  )
}))

p7_cond_dens <- ggplot(dens_df, aes(x = x, y = density, color = Age_Group)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "",
    x = "Survival Time (Months)",
    y = "Estimated Conditional Density f(x | z)",
    color = ""
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Display all plots
print(p1_scatter)
print(p2_theta_curve)
print(p3_cond_mean)
print(p4_cv_curve)
print(p5_sens_bw)
print(p6_boundary)
print(p7_cond_dens)

# Plot 9: Histogram of Survival Time (X) with Local Constant (p=0),
# Local Linear (p=1), Local Quadratic (p=2), and Parametric Densities
# (Log-Linear Lindley omitted)
x_grid <- seq(min(df$X), max(df$X), length.out = 300)
z_mean <- mean(df$Z)

# 1. Local Constant Lindley (p = 0) at z_mean
th_loc0 <- fit_local_lindley_point(df$X, df$Z, z0 = z_mean, h = h_opt_0, p = 0)$theta0
dens_loc0 <- exp(lindley_logpdf(x_grid, th_loc0))

# 2. Local Linear Lindley (p = 1) at z_mean
th_loc1 <- fit_local_lindley_point(df$X, df$Z, z0 = z_mean, h = h_opt_1, p = 1)$theta0
dens_loc1 <- exp(lindley_logpdf(x_grid, th_loc1))

# 3. Local Quadratic Lindley (p = 2) at z_mean
th_loc2 <- fit_local_lindley_point(df$X, df$Z, z0 = z_mean, h = h_opt_2, p = 2)$theta0
dens_loc2 <- exp(lindley_logpdf(x_grid, th_loc2))

# 4. Constant Lindley (marginal fit)
dens_const <- exp(lindley_logpdf(x_grid, opt_const$minimum))

# 5. Exponential Regression at z_mean
mu_exp_mean <- as.numeric(predict(fit_exp, newdata = data.frame(Z = z_mean), type = "response"))
dens_exp <- dexp(x_grid, rate = 1 / mu_exp_mean)

# 6. Gamma Regression at z_mean
mu_gam_mean <- as.numeric(predict(fit_gam, newdata = data.frame(Z = z_mean), type = "response"))
shape_gam <- 1 / summary(fit_gam)$dispersion
dens_gam <- dgamma(x_grid, shape = shape_gam, scale = mu_gam_mean / shape_gam)

# 7. Weibull Regression at z_mean
sc_wei_mean <- as.numeric(predict(fit_wei, newdata = data.frame(Z = z_mean), type = "response"))
shape_wei <- 1 / fit_wei$scale
dens_wei <- dweibull(x_grid, shape = shape_wei, scale = sc_wei_mean)

hist_dens_df <- rbind(
  data.frame(x = x_grid, density = dens_loc0,  Model = "Local Constant Lindley (p = 0)"),
  data.frame(x = x_grid, density = dens_loc1,  Model = "Local Linear Lindley (p = 1)"),
  data.frame(x = x_grid, density = dens_loc2,  Model = "Local Quadratic Lindley (p = 2)"),
  data.frame(x = x_grid, density = dens_const, Model = "Constant Lindley"),
  data.frame(x = x_grid, density = dens_exp,   Model = "Exponential Regression"),
  data.frame(x = x_grid, density = dens_gam,   Model = "Gamma Regression"),
  data.frame(x = x_grid, density = dens_wei,   Model = "Weibull Regression")
)

p9_hist_models <- ggplot() +
  geom_histogram(
    data = df,
    aes(x = X, y = after_stat(density)),
    bins = 20,
    fill = "gray85",
    color = "gray40",
    linewidth = 0.4
  ) +
  geom_line(
    data = hist_dens_df,
    aes(x = x, y = density, color = Model, linetype = Model),
    linewidth = 1.05
  ) +
  labs(
    title = "",
    x = "Survival Time (Months)",
    y = "Density",
    color = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p9_hist_models)

