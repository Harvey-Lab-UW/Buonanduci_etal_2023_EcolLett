# Scaling analysis
# Evaluate whether relationship btwn fire size and spatial metrics
# varies by region

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
  mutate(Region = factor(Region)) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

# Pull out unique regions/regimes
fire_metrics_high <- filter(fire_metrics, Fire_Regime == "High")
fire_metrics_mix <- filter(fire_metrics, Fire_Regime == "Mixed")
fire_metrics_low <- filter(fire_metrics, Fire_Regime == "Low")

# Create facet labels
FRG_names <- c(
  `High` = "Regime: High",
  `Mixed` = "Regime: Mixed",
  `Low` = "Regime: Low"
)


# Define parameters and functions -------------

# Taus (quantiles) to evaluate
taus <- c(0.05,0.5,0.95)

# Confidence intervals
conf <- 0.95
z <- -qnorm((1-conf)/2)

# Function to create data frame of predictor variables for region-specific curves
new_data_fn <- function(df){
  t1 <- tibble(log_fire_area = seq(min(df$log_fire_area[df$Region == "Northern Rockies"]), 
                                      max(df$log_fire_area[df$Region == "Northern Rockies"]), by = 0.1),
               Region = "Northern Rockies")
  t2 <- tibble(log_fire_area = seq(min(df$log_fire_area[df$Region == "Pacific Northwest"]), 
                                      max(df$log_fire_area[df$Region == "Pacific Northwest"]), by = 0.1),
               Region = "Pacific Northwest")
  t <- bind_rows(t1, t2) %>%
    mutate(Region = factor(Region))
  return(t)
} 

# Function to make predictions for region-specific curves
predict_region_fn <- function(models, df){
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
  bind_rows(pred_list) %>% mutate(Fire_Regime = df$Fire_Regime[1])
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
new_data_dif_fn <- function(df, var){
  # Find overlapping range of fire sizes for both regions
  size_range <- df %>%
    drop_na(.data[[var]]) %>%
    group_by(Region) %>%
    dplyr::summarize(fire_min = min(log_fire_area), fire_max = max(log_fire_area))
  # Create new_data
  new_data <- expand_grid(log_fire_area = seq(max(size_range$fire_min),
                                                 min(size_range$fire_max), by = 0.1),
                          Region = c("Northern Rockies", "Pacific Northwest"))
  return(new_data)
}

# Function to make predictions for difference between curves
predict_region_diff_fn <- function(models, df, var){
  # List to store predictions
  dif_list <- list()
  # Dataframe of predictor variables
  new_df <- new_data_dif_fn(df, var)
  # Loop through tau-specific models and make predictions
  for (i in 1:length(taus)){
    dif_list[[i]] <- smooth_diff(models[[i]], newdata = new_df,
                                 f1 = "Pacific Northwest", f2 = "Northern Rockies", var = "Region", alpha = 1-conf) %>%
      mutate(tau = taus[i]) %>%
      mutate(log_fire_area = unique(new_df$log_fire_area))
  }
  bind_rows(dif_list) %>% mutate(Fire_Regime = df$Fire_Regime[1])
}



# Area-wtd mean patch size ----

# Fit quantile- and regime-specific models
mod_high <- list()
mod_mix <- list()
mod_low <- list()

for(i in 1:length(taus)){
  mod_high[[i]] <- gcrq(log_patch_area_AW_mean ~ Region + 
                          ps(log_fire_area, shared.pen = TRUE, monotone = 1, by = Region),
                        lambda0 = 5, tau = taus[i], data = fire_metrics_high)
  mod_mix[[i]] <- gcrq(log_patch_area_AW_mean ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, monotone = 1, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_mix)
  mod_low[[i]] <- gcrq(log_patch_area_AW_mean ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, monotone = 1, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_low)
}

# Make predictions
pred_high <- predict_region_fn(mod_high, fire_metrics_high)
pred_mix <- predict_region_fn(mod_mix, fire_metrics_mix)
pred_low <- predict_region_fn(mod_low, fire_metrics_low)

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Region, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

# Plotting
p1 = ggplot( filter(pred_wide, Fire_Regime == "Low" | Fire_Regime == "Mixed") ) +
  facet_grid(~ Fire_Regime, labeller = as_labeller(FRG_names)) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Patch size:\nArea-wtd mean (ha)") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  scale_color_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  scale_fill_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  geom_point(data =  filter(fire_metrics, Fire_Regime == "Low" | Fire_Regime == "Mixed") , 
             aes(x = log_fire_area, y = log_patch_area_AW_mean, color = Region), size = 0.6, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Region), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Region)) +
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
dif_high <- predict_region_diff_fn(mod_high, fire_metrics_high, "log_patch_area_AW_mean")
dif_mix <- predict_region_diff_fn(mod_mix, fire_metrics_mix, "log_patch_area_AW_mean")
dif_low <- predict_region_diff_fn(mod_low, fire_metrics_low, "log_patch_area_AW_mean")

dif <- bind_rows(dif_high, dif_mix, dif_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05))) %>%
  filter(Fire_Regime == "Low" | Fire_Regime == "Mixed")

# Plotting
p2 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ Fire_Regime, labeller = as_labeller(FRG_names)) +
  labs(y = "Difference*") +
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
        strip.text = element_text(size = 9),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")


# Patch size distribution: beta ----

# Fit quantile- and regime-specific models
mod_high <- list()
mod_mix <- list()
mod_low <- list()

for(i in 1:length(taus)){
  mod_high[[i]] <- gcrq(beta ~ Region + 
                          ps(log_fire_area, shared.pen = TRUE, by = Region),
                        lambda0 = 5, tau = taus[i], data = fire_metrics_high)
  mod_mix[[i]] <- gcrq(beta ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_mix)
  mod_low[[i]] <- gcrq(beta ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_low)
}

# Make predictions
pred_high <- predict_region_fn(mod_high, fire_metrics_high)
pred_mix <- predict_region_fn(mod_mix, fire_metrics_mix)
pred_low <- predict_region_fn(mod_low, fire_metrics_low)

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Region, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

# Plotting
p3 = ggplot( filter(pred_wide, Fire_Regime == "Low" | Fire_Regime == "Mixed") ) +
  facet_grid(~ Fire_Regime) +
  scale_y_continuous(name = "Patch size:\n\u03B2 parameter") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  scale_color_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  scale_fill_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  geom_point(data =  filter(fire_metrics, Fire_Regime == "Low" | Fire_Regime == "Mixed") , 
             aes(x = log_fire_area, y = beta, color = Region), size = 0.6, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Region), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Region)) +
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

# Predictions for differences between curves
dif_high <- predict_region_diff_fn(mod_high, fire_metrics_high, "beta")
dif_mix <- predict_region_diff_fn(mod_mix, fire_metrics_mix, "beta")
dif_low <- predict_region_diff_fn(mod_low, fire_metrics_low, "beta")

dif <- bind_rows(dif_high, dif_mix, dif_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05))) %>%
  filter(Fire_Regime == "Low" | Fire_Regime == "Mixed")

# Plotting
p4 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ Fire_Regime) +
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

# Fit quantile- and regime-specific models
mod_high <- list()
mod_mix <- list()
mod_low <- list()

for(i in 1:length(taus)){
  mod_high[[i]] <- gcrq(psi ~ Region + 
                          ps(log_fire_area, shared.pen = TRUE, by = Region),
                        lambda0 = 5, tau = taus[i], data = fire_metrics_high)
  mod_mix[[i]] <- gcrq(psi ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_mix)
  mod_low[[i]] <- gcrq(psi ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_low)
}

# Make predictions
pred_high <- predict_region_fn(mod_high, fire_metrics_high)
pred_mix <- predict_region_fn(mod_mix, fire_metrics_mix)
pred_low <- predict_region_fn(mod_low, fire_metrics_low)

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Region, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

# Plotting
p5 = ggplot( filter(pred_wide, Fire_Regime == "Low" | Fire_Regime == "Mixed") ) +
  facet_grid(~ Fire_Regime) +
  scale_y_continuous(limits = c(-0.21, 1),
                     name = "Patch size:\n\u03C8 parameter") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  scale_color_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  scale_fill_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  geom_point(data =  filter(fire_metrics, Fire_Regime == "Low" | Fire_Regime == "Mixed") , 
             aes(x = log_fire_area, y = psi, color = Region), size = 0.6, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Region), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Region)) +
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

# Predictions for differences between curves
dif_high <- predict_region_diff_fn(mod_high, fire_metrics_high, "psi")
dif_mix <- predict_region_diff_fn(mod_mix, fire_metrics_mix, "psi")
dif_low <- predict_region_diff_fn(mod_low, fire_metrics_low, "psi")

dif <- bind_rows(dif_high, dif_mix, dif_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05))) %>%
  filter(Fire_Regime == "Low" | Fire_Regime == "Mixed")

# Plotting
p6 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ Fire_Regime) +
  labs(y = "Difference") +
  scale_color_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_fill_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  coord_cartesian(ylim = c(-0.4,0.4)) +
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

# Fit quantile- and regime-specific models
mod_high <- list()
mod_mix <- list()
mod_low <- list()

for(i in 1:length(taus)){
  mod_high[[i]] <- gcrq(log_total_core ~ Region + 
                          ps(log_fire_area, shared.pen = TRUE, monotone = 1, by = Region),
                        lambda0 = 5, tau = taus[i], data = fire_metrics_high)
  mod_mix[[i]] <- gcrq(log_total_core ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, monotone = 1, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_mix)
  mod_low[[i]] <- gcrq(log_total_core ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, monotone = 1, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_low)
}

# Make predictions
pred_high <- predict_region_fn(mod_high, fire_metrics_high)
pred_mix <- predict_region_fn(mod_mix, fire_metrics_mix)
pred_low <- predict_region_fn(mod_low, fire_metrics_low)

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Region, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

# Plotting
p7 = ggplot( filter(pred_wide, Fire_Regime == "Low" | Fire_Regime == "Mixed") ) +
  facet_grid(~ Fire_Regime) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Patch structure:\nTotal core area (ha)") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  scale_color_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  scale_fill_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  geom_point(data =  filter(fire_metrics, Fire_Regime == "Low" | Fire_Regime == "Mixed") , 
             aes(x = log_fire_area, y = log_total_core, color = Region), size = 0.6, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Region), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Region)) +
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

# Predictions for differences between curves
dif_high <- predict_region_diff_fn(mod_high, fire_metrics_high, "log_total_core")
dif_mix <- predict_region_diff_fn(mod_mix, fire_metrics_mix, "log_total_core")
dif_low <- predict_region_diff_fn(mod_low, fire_metrics_low, "log_total_core")

dif <- bind_rows(dif_high, dif_mix, dif_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05))) %>%
  filter(Fire_Regime == "Low" | Fire_Regime == "Mixed")

# Plotting
p8 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ Fire_Regime) +
  labs(y = "Difference*") +
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




# Seed decay coefficient ----

# Fit quantile- and regime-specific models

# In high-severity regime,
# issue with singularity in matrix for quantile regression @ quantile 0.95
# adjust slightly for model fitting
taus_hs <- c(0.05,0.5,0.949)

mod_high <- list()
mod_mix <- list()
mod_low <- list()

for(i in 1:length(taus)){
  mod_high[[i]] <- gcrq(log_SDC ~ Region + 
                          ps(log_fire_area, shared.pen = TRUE, monotone = -1, by = Region),
                        lambda0 = 5, tau = taus_hs[i], data = fire_metrics_high)
  mod_mix[[i]] <- gcrq(log_SDC ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, monotone = -1, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_mix)
  mod_low[[i]] <- gcrq(log_SDC ~ Region +
                         ps(log_fire_area, shared.pen = TRUE, monotone = -1, by = Region),
                       lambda0 = 5, tau = taus[i], data = fire_metrics_low)
}

# Make predictions
pred_high <- predict_region_fn(mod_high, fire_metrics_high) %>%
  mutate(tau = round(tau, 2))
pred_mix <- predict_region_fn(mod_mix, fire_metrics_mix)
pred_low <- predict_region_fn(mod_low, fire_metrics_low)

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

pred_wide <- pred %>%
  select(tau, fit, log_fire_area, Region, Fire_Regime) %>%
  pivot_wider(names_from = tau, values_from = fit)

# Plotting
p9 = ggplot( filter(pred_wide, Fire_Regime == "Low" | Fire_Regime == "Mixed") ) +
  facet_grid(~ Fire_Regime) +
  scale_y_continuous(breaks = c(-3, -2.5, -2, -1.5),
                     labels = c("0.001", "0.003", "0.01", "0.03"),
                     name = "Patch structure:\nSDC parameter") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  scale_color_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  scale_fill_manual(values = c("tomato3", "gray20"), labels = c("Northern\nRockies", "Pacific\nNorthwest")) +
  geom_point(data =  filter(fire_metrics, Fire_Regime == "Low" | Fire_Regime == "Mixed") , 
             aes(x = log_fire_area, y = log_SDC, color = Region), size = 0.6, alpha = 0.2) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Region), alpha = 0.15) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Region), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Region)) +
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
dif_high <- predict_region_diff_fn(mod_high, fire_metrics_high, "log_SDC")
dif_mix <- predict_region_diff_fn(mod_mix, fire_metrics_mix, "log_SDC")
dif_low <- predict_region_diff_fn(mod_low, fire_metrics_low, "log_SDC")

dif <- bind_rows(dif_high, dif_mix, dif_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau, levels = c(0.95, 0.5, 0.05))) %>%
  filter(Fire_Regime == "Low" | Fire_Regime == "Mixed")

# Plotting
p10 = ggplot(dif, aes(x = log_fire_area, y = diff, color = tau, fill = tau)) +
  facet_grid(~ Fire_Regime) +
  labs(y = "Difference*") +
  scale_color_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_fill_viridis_d(name = "Quantile", option = "mako", end = 0.7) +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  coord_cartesian(ylim = c(-0.4,0.4)) +
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



# Combine all plots ---------

g1 <- plot_grid(p1, p3, p5, p7, p9,
          nrow = 5, align = "v", rel_heights = c(1.2, 1, 1, 1, 1.7))

g2 <- plot_grid(p2, p4, p6, p8, p10,
                nrow = 5, align = "v", rel_heights = c(1.2, 1, 1, 1, 1.7))

plot_grid(g1, g2, nrow = 1, rel_widths = c(1.2, 1), 
          labels = c("(a)", "(b)"), label_size = 10)


# Export figure
ggsave(paste0("Figures/figureS2_scaling_regions.png"),
       width = 173, height = 165, units = "mm", dpi = 600)








