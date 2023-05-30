# Scaling analysis
# Evaluate whether relationship btwn fire size and spatial metrics
# varies by fire regime

# Load packages
require(tidyverse)
require(quantreg)
require(quantregGrowth)
require(cowplot)

# Load functions
# Function to extract linear predictor matrix from gcrq model object
source("Scripts/predict.gcrq.lpmatrix.R")


# Load data--------------

# Load landscape metrics data file
fire_metrics <- read_csv("Data/fire_metrics.csv")

# Create factors
fire_metrics <- fire_metrics %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))



# Define parameters and functions -------------

# Taus (quantiles) to evaluate
taus <- c(0.05, 0.5, 0.95)

# Confidence intervals
conf <- 0.95
z <- -qnorm((1-conf)/2)

# Function to create data frame of predictor variables for regime-specific curves
new_data_fn <- function(df){
  t1 <- tibble(log_fire_area = seq(min(df$log_fire_area[df$Fire_Regime == "Low"]), 
                                      max(df$log_fire_area[df$Fire_Regime == "Low"]), by = 0.1),
               Fire_Regime = "Low")
  t2 <- tibble(log_fire_area = seq(min(df$log_fire_area[df$Fire_Regime == "Mixed"]), 
                                      max(df$log_fire_area[df$Fire_Regime == "Mixed"]), by = 0.1),
               Fire_Regime = "Mixed")
  t3 <- tibble(log_fire_area = seq(min(df$log_fire_area[df$Fire_Regime == "High"]), 
                                      max(df$log_fire_area[df$Fire_Regime == "High"]), by = 0.1),
               Fire_Regime = "High")
  t <- bind_rows(t1, t2, t3) %>%
    mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))
  return(t)
} 

# Function to make predictions for regime-specific curves
predict_regime_fn <- function(models, df){
  # List to store predictions
  pred_list <- list()
  # Dataframe of predictor variables
  new_df <- new_data_fn(df)
  # Loop through tau-specific models and make predictions
  for (i in 1:length(taus)){
    pred <- predict.gcrq(models[[i]], newdata = as.data.frame(new_df), se.fit = TRUE)
    pred_list[[i]] <- tibble(tau = taus[i], fit = pred$fit, se.fit = pred$se.fit) %>%
      mutate(lwr = fit-z*se.fit) %>% mutate(upr = fit+z*se.fit) %>%
      bind_cols(new_df)
  }
  bind_rows(pred_list)
}

# Functions to evaluate difference between smooth curves in each model:

# Function to calculate difference between smooth regression curves
# Adapted from https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/
smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05) {
  xp <- predict.gcrq.lpmatrix(model, newdata = as.data.frame(newdata))
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other factors
  X[, ! (c1 | c2)] <- 0

  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model)[[1]]) * X))
  crit <- qnorm(alpha/2, lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = as.numeric(dif),
             se = as.numeric(se),
             upper = as.numeric(upr),
             lower = as.numeric(lwr))
}

# Function to create new data frame for *difference* predictions
new_data_dif_fn <- function(df, var, f1, f2){
  # Find overlapping range of fire sizes for both regimes
  size_range <- df %>%
    drop_na(.data[[var]]) %>%
    filter(Fire_Regime == f1 | Fire_Regime == f2) %>%
    group_by(Fire_Regime) %>%
    dplyr::summarize(fire_min = min(log_fire_area), fire_max = max(log_fire_area))
  # Create new_data
  new_data <- expand_grid(log_fire_area = seq(max(size_range$fire_min),
                                                 min(size_range$fire_max), by = 0.1),
                          Fire_Regime = c("Low", "Mixed", "High"))
  return(new_data)
}

# Function to make predictions for difference between curves
predict_regime_diff_fn <- function(models, df, var, f1, f2){
  # List to store predictions
  dif_list <- list()
  # Dataframe of predictor variables
  new_df <- new_data_dif_fn(df, var, f1, f2)
  # Loop through tau-specific models and make predictions
  for (i in 1:length(taus)){
    dif_list[[i]] <- smooth_diff(models[[i]], newdata = new_df,
                                 f1 = f1, f2 = f2, var = "Fire_Regime", alpha = 1-conf) %>%
      mutate(tau = taus[i]) %>%
      mutate(log_fire_area = unique(new_df$log_fire_area))
  }
  bind_rows(dif_list)
}




# Area-wtd mean patch size ----

# Fit quantile-specific models
mod <- list()

for(i in 1:length(taus)){
  mod[[i]] <- gcrq(log_patch_area_AW_mean ~ Fire_Regime + 
                     ps(log_fire_area, monotone = 1, shared.pen = TRUE, by = Fire_Regime),
                   lambda0 = 5, tau = taus[i], data = fire_metrics)
}

# Make predictions
pred <- predict_regime_fn(mod, fire_metrics)

# Plotting
pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

p1 = ggplot(pred_wide) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Patch size:\nArea-wtd mean (ha)") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  geom_point(data = fire_metrics, 
             aes(x = log_fire_area, y = log_patch_area_AW_mean, color = Fire_Regime), size = 0.6, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  ggtitle(" ") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 9),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")


# Make predictions for differences between curves
dif_high_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "log_patch_area_AW_mean", "Low", "High")
dif_high_mix <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "log_patch_area_AW_mean", "Mixed", "High")
dif_mix_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                      "log_patch_area_AW_mean", "Low", "Mixed")
dif <- bind_rows(dif_high_low, dif_high_mix, dif_mix_low) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05)))


# Plotting 
p2 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ pair) +
  labs(y = "Difference*") +
  scale_color_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_fill_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  coord_cartesian(ylim = c(-1.5,1)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 9),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")



# Patch size distribution: beta ----

# Fit quantile-specific models
mod <- list()

for(i in 1:length(taus)){
  mod[[i]] <- gcrq(beta ~ Fire_Regime + 
                     ps(log_fire_area, shared.pen = TRUE, by = Fire_Regime),
                   lambda0 = 5, tau = taus[i], data = fire_metrics)
}

# Make predictions
pred <- predict_regime_fn(mod, fire_metrics)

# Plotting 
pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

p3 = ggplot(pred_wide) +
  scale_y_continuous(name = "Patch size:\n\u03B2 parameter") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  geom_point(data = fire_metrics, 
             aes(x = log_fire_area, y = beta, color = Fire_Regime), size = 0.6, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 9),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")


# Predictions for differences between curves
dif_high_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "beta", "Low", "High")
dif_high_mix <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "beta", "Mixed", "High")
dif_mix_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                      "beta", "Low", "Mixed")
dif <- bind_rows(dif_high_low, dif_high_mix, dif_mix_low) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05)))


# Plotting 
p4 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ pair) +
  labs(y = "Difference") +
  scale_color_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_fill_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  coord_cartesian(ylim = c(-1.5,1.5)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")



# Patch size distribution: psi ----

# Fit quantile-specific models
mod <- list()

for(i in 1:length(taus)){
  mod[[i]] <- gcrq(psi ~ Fire_Regime + 
                     ps(log_fire_area, shared.pen = TRUE, by = Fire_Regime),
                   lambda0 = 5, tau = taus[i], data = fire_metrics)
}

# Make predictions
pred <- predict_regime_fn(mod, fire_metrics)

# Plotting 
pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

p5 = ggplot(pred_wide) +
  scale_y_continuous(limits = c(-0.21, 1),
                     name = "Patch size:\n\u03C8 parameter") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  geom_point(data = fire_metrics, 
             aes(x = log_fire_area, y = psi, color = Fire_Regime), size = 0.6, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 9),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")


# Predictions for differences between curves
dif_high_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "psi", "Low", "High")
dif_high_mix <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "psi", "Mixed", "High")
dif_mix_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                      "psi", "Low", "Mixed")
dif <- bind_rows(dif_high_low, dif_high_mix, dif_mix_low) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05)))


# Plotting 
p6 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ pair) +
  labs(y = "Difference") +
  scale_color_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_fill_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")



# Core area ----

# Fit quantile-specific models
mod <- list()

for(i in 1:length(taus)){
  mod[[i]] <- gcrq(log_total_core ~ Fire_Regime + 
                     ps(log_fire_area, monotone = 1, shared.pen = TRUE, by = Fire_Regime),
                   lambda0 = 5, tau = taus[i], data = fire_metrics)
}

# Make predictions
pred <- predict_regime_fn(mod, fire_metrics)

# Plotting
pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

p7 = ggplot(pred_wide) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Patch structure:\nTotal core area (ha)") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  geom_point(data = fire_metrics, 
             aes(x = log_fire_area, y = log_total_core, color = Fire_Regime), size = 0.6, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 9),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")


# Predictions for differences between curves
dif_high_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "log_total_core", "Low", "High")
dif_high_mix <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "log_total_core", "Mixed", "High")
dif_mix_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                      "log_total_core", "Low", "Mixed")
dif <- bind_rows(dif_high_low, dif_high_mix, dif_mix_low) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05)))


# Plotting 
p8 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ pair) +
  labs(y = "Difference*") +
  scale_color_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_fill_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  coord_cartesian(ylim = c(-3,1)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")



# Seed decay coefficient ----

# Fit quantile-specific models
mod <- list()

for(i in 1:length(taus)){
  mod[[i]] <- gcrq(log_SDC ~ Fire_Regime + 
                     ps(log_fire_area, monotone = -1, shared.pen = TRUE, by = Fire_Regime),
                   lambda0 = 5, tau = taus[i], data = fire_metrics)
}

# Make predictions
pred <- predict_regime_fn(mod, fire_metrics)

# Plotting
pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

p9 = ggplot(pred_wide) +
  scale_y_continuous(breaks = c(-3, -2.5, -2, -1.5),
                     labels = c("0.001", "0.003", "0.01", "0.03"),
                     name = "Patch structure:\nSDC parameter") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  scale_fill_viridis_d(name = NULL,
                       direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(name = NULL,
                        direction = -1, option = "inferno", end = 0.8) +
  geom_point(data = fire_metrics, 
             aes(x = log_fire_area, y = log_SDC, color = Fire_Regime), size = 0.6, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = "bottom")


# Predictions for differences between curves
dif_high_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "log_SDC", "Low", "High")
dif_high_mix <- predict_regime_diff_fn(mod, fire_metrics, 
                                       "log_SDC", "Mixed", "High")
dif_mix_low <- predict_regime_diff_fn(mod, fire_metrics, 
                                      "log_SDC", "Low", "Mixed")
dif <- bind_rows(dif_high_low, dif_high_mix, dif_mix_low) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05)))


# Plotting 
p10 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ pair) +
  labs(y = "Difference*") +
  scale_color_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_fill_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = "bottom")



# Combine all plots

g1 <- plot_grid(p1, p3, p5, p7, p9,
                nrow = 5, align = "v", rel_heights = c(1.2, 1, 1, 1, 1.7))

g2 <- plot_grid(p2, p4, p6, p8, p10,
                nrow = 5, align = "v", rel_heights = c(1.2, 1, 1, 1, 1.7))

plot_grid(g1, g2, nrow = 1, rel_widths = c(1.2, 2), 
          labels = c("(a)", "(b)"), label_size = 10)


# Export figure
ggsave("Figures/figureS1_scaling_regimes.png", bg = "white",
       width = 173, height = 165, units = "mm", dpi = 600)




