# Scaling analysis
# Figure illustrating marginal effect of Year on scaling relationships

# Load packages
require(tidyverse)
require(quantreg)
require(quantregGrowth)
require(cowplot)

# Load data --------------

# Load landscape metrics data file
fire_metrics <- read_csv("Data/fire_metrics.csv")

# Create factors
fire_metrics <- fire_metrics %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Mean-center and standardize year for use as linear predictor
fire_metrics <- fire_metrics %>%
  mutate(Year_std = (Year - mean(Year))/sd(Year))

# Pull out unique fire regimes
fire_metrics_high <- filter(fire_metrics, Fire_Regime == "High")
fire_metrics_mix <- filter(fire_metrics, Fire_Regime == "Mixed")
fire_metrics_low <- filter(fire_metrics, Fire_Regime == "Low")

# Create facet labels
FRG_names <- c(
  `High` = "Regime: High",
  `Mixed` = "Regime: Mixed",
  `Low` = "Regime: Low"
)

# Prediction function
predict_year_fn <- function(model, df){
  # Create data frame of fire sizes and years
  new_df <- cross_df(list(log_fire_area = seq(2.5, max(df$log_fire_area), by = 0.1),
                          Year = c(1985:2020)))
  # Standardize Year
  new_df <- new_df %>%
    mutate(Year_std = (Year - mean(fire_metrics$Year)) / sd(fire_metrics$Year))
  
  # Make predictions
  pred <- predict.gcrq(model, newdata = as.data.frame(new_df), se.fit = TRUE) %>%
    as_tibble() %>% bind_cols(new_df) %>% 
    mutate(Region = df$Region[1]) %>% mutate(Fire_Regime = df$Fire_Regime[1])
  
  # If multiple taus, unpack fit and se columns
  if ( !is.null(ncol(pred$fit)) ){
    pred <- pred %>%
      mutate(fit = as_tibble(fit)) %>% mutate(se.fit = as_tibble(se.fit)) %>%
      unpack(cols = fit, names_sep = "_") %>%
      unpack(cols = se.fit, names_sep = "_")
  }
}

# Function to extract p-values for linear year terms
extract_p <- function(model, linterm){
  n.tau <- ncol(as.matrix(model$coefficients))
  list.vcov <- vcov.gcrq(model)
  est <- c()
  se <- c()
  p_val <- c()
  for (j in 1:n.tau) {
    est <- c(est, as.matrix(model$coefficients)[linterm, j])
    se <- c(se, sqrt(diag(list.vcov[[j]]))[linterm])
    p_val <- c(p_val, pchisq((est[j]/se[j])^2, df = 1, lower.tail = FALSE))
  }
  return(tibble(tau = model$taus, 
                est = est, se = se, p_val = p_val))
}


# Area-wtd mean patch size ----

# Fit quantile regression models: additive smooth term for year
mod_smoo_high <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + ps(Year),
                      lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_smoo_mix <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_smoo_low <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Fit quantile regression models: additive liear term for year
mod_lin_high <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + Year_std,
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_lin_mix <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_lin_low <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions for additive smooth
pred_smoo_high <- predict_year_fn(mod_smoo_high, fire_metrics_high)
pred_smoo_mix <- predict_year_fn(mod_smoo_mix, fire_metrics_mix)
pred_smoo_low <- predict_year_fn(mod_smoo_low, fire_metrics_low)

pred_smoo <- bind_rows(pred_smoo_high, pred_smoo_mix, pred_smoo_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Predictions for additive linear (with p-values)
pred_lin_high <- predict_year_fn(mod_lin_high, fire_metrics_high) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_high, "Year_std") ) # join p-values for linear terms
pred_lin_mix <- predict_year_fn(mod_lin_mix, fire_metrics_mix) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_mix, "Year_std") ) # join p-values for linear terms
pred_lin_low <- predict_year_fn(mod_lin_low, fire_metrics_low) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_low, "Year_std") ) # join p-values for linear terms

pred_lin <- bind_rows(pred_lin_high, pred_lin_mix, pred_lin_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))


# Plotting
p1 <- ggplot( filter(pred_smoo, log_fire_area == 3.5) ) +
  facet_grid( ~ Fire_Regime, labeller = as_labeller(FRG_names)) +
  scale_y_continuous(limits = range(fire_metrics$log_patch_area_AW_mean, na.rm = T), 
                     breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Area-weighted\nmean (ha)") +
  scale_x_continuous(name = "Year") +
  geom_ribbon(aes(x = Year, ymin = fit_0.05, ymax = fit_0.95, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = Year, ymin = fit_0.25, ymax = fit_0.75, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = Year, y = fit_0.05, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.25, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.75, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.95, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.5, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val >= 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.2, lty = 3) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val < 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.5, lty = 3) +
  ggtitle("Patch size metrics") +
  theme_bw() +
  theme(plot.title = element_text(size = 9, vjust = 3, face = "italic"),
        plot.title.position = "plot",
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 9),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")


# Patch distribution: beta ----

# Fit quantile regression models: additive smooth term for year
mod_smoo_high <- gcrq(beta ~ ps(log_fire_area) + ps(Year),
                      lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_smoo_mix <- gcrq(beta ~ ps(log_fire_area) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_smoo_low <- gcrq(beta ~ ps(log_fire_area) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Fit quantile regression models: additive linear term for year
mod_lin_high <- gcrq(beta ~ ps(log_fire_area) + Year_std,
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_lin_mix <- gcrq(beta ~ ps(log_fire_area) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_lin_low <- gcrq(beta ~ ps(log_fire_area) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions for additive smooth
pred_smoo_high <- predict_year_fn(mod_smoo_high, fire_metrics_high)
pred_smoo_mix <- predict_year_fn(mod_smoo_mix, fire_metrics_mix)
pred_smoo_low <- predict_year_fn(mod_smoo_low, fire_metrics_low)

pred_smoo <- bind_rows(pred_smoo_high, pred_smoo_mix, pred_smoo_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Predictions for additive linear (with p-values)
pred_lin_high <- predict_year_fn(mod_lin_high, fire_metrics_high) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_high, "Year_std") ) # join p-values for linear terms
pred_lin_mix <- predict_year_fn(mod_lin_mix, fire_metrics_mix) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_mix, "Year_std") ) # join p-values for linear terms
pred_lin_low <- predict_year_fn(mod_lin_low, fire_metrics_low) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_low, "Year_std") ) # join p-values for linear terms

pred_lin <- bind_rows(pred_lin_high, pred_lin_mix, pred_lin_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))


# Plotting
p2 <- ggplot( filter(pred_smoo, log_fire_area == 3.5) ) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(limits = range(fire_metrics$beta, na.rm = T), 
                     name = "\u03B2\nparameter") +
  scale_x_continuous(name = "Year") +
  geom_ribbon(aes(x = Year, ymin = fit_0.05, ymax = fit_0.95, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = Year, ymin = fit_0.25, ymax = fit_0.75, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = Year, y = fit_0.05, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.25, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.75, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.95, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.5, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val >= 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.2, lty = 3) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val < 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.5, lty = 3) +
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



# Patch distribution: psi ----

# Fit quantile regression models: additive smooth term for year
mod_smoo_high <- gcrq(psi ~ ps(log_fire_area) + ps(Year),
                      lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_smoo_mix <- gcrq(psi ~ ps(log_fire_area) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_smoo_low <- gcrq(psi ~ ps(log_fire_area) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Fit quantile regression models: additive linear term for year
mod_lin_high <- gcrq(psi ~ ps(log_fire_area) + Year_std,
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_lin_mix <- gcrq(psi ~ ps(log_fire_area) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_lin_low <- gcrq(psi ~ ps(log_fire_area) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions for additive smooth
pred_smoo_high <- predict_year_fn(mod_smoo_high, fire_metrics_high)
pred_smoo_mix <- predict_year_fn(mod_smoo_mix, fire_metrics_mix)
pred_smoo_low <- predict_year_fn(mod_smoo_low, fire_metrics_low)

pred_smoo <- bind_rows(pred_smoo_high, pred_smoo_mix, pred_smoo_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Predictions for additive linear (with p-values)
pred_lin_high <- predict_year_fn(mod_lin_high, fire_metrics_high) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_high, "Year_std") ) # join p-values for linear terms
pred_lin_mix <- predict_year_fn(mod_lin_mix, fire_metrics_mix) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_mix, "Year_std") ) # join p-values for linear terms
pred_lin_low <- predict_year_fn(mod_lin_low, fire_metrics_low) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_low, "Year_std") ) # join p-values for linear terms

pred_lin <- bind_rows(pred_lin_high, pred_lin_mix, pred_lin_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))


# Plotting
p3 <- ggplot( filter(pred_smoo, log_fire_area == 3.5) ) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(limits = c(-0.21, 1),
                     name = "\u03C8\nparameter") +
  scale_x_continuous(name = "Year") +
  geom_ribbon(aes(x = Year, ymin = fit_0.05, ymax = fit_0.95, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = Year, ymin = fit_0.25, ymax = fit_0.75, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = Year, y = fit_0.05, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.25, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.75, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.95, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.5, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val >= 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.2, lty = 3) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val < 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.5, lty = 3) +
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

# Fit quantile regression models: additive smooth term for year
mod_smoo_high <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1) + ps(Year),
                      lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_smoo_mix <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_smoo_low <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Fit quantile regression models: additive linear term for year
mod_lin_high <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1) + Year_std,
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_lin_mix <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_lin_low <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions for additive smooth
pred_smoo_high <- predict_year_fn(mod_smoo_high, fire_metrics_high)
pred_smoo_mix <- predict_year_fn(mod_smoo_mix, fire_metrics_mix)
pred_smoo_low <- predict_year_fn(mod_smoo_low, fire_metrics_low)

pred_smoo <- bind_rows(pred_smoo_high, pred_smoo_mix, pred_smoo_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Predictions for additive linear (with p-values)
pred_lin_high <- predict_year_fn(mod_lin_high, fire_metrics_high) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_high, "Year_std") ) # join p-values for linear terms
pred_lin_mix <- predict_year_fn(mod_lin_mix, fire_metrics_mix) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_mix, "Year_std") ) # join p-values for linear terms
pred_lin_low <- predict_year_fn(mod_lin_low, fire_metrics_low) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_low, "Year_std") ) # join p-values for linear terms

pred_lin <- bind_rows(pred_lin_high, pred_lin_mix, pred_lin_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))


# Plotting
p4 <- ggplot( filter(pred_smoo, log_fire_area == 3.5) ) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(limits = range(fire_metrics$log_total_core, na.rm = T), 
                     breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Total core\narea (ha)") +
  scale_x_continuous(name = "Year") +
  geom_ribbon(aes(x = Year, ymin = fit_0.05, ymax = fit_0.95, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = Year, ymin = fit_0.25, ymax = fit_0.75, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = Year, y = fit_0.05, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.25, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.75, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.95, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.5, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val >= 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.2, lty = 3) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val < 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.5, lty = 3) +
  ggtitle("Patch structure metrics") +
  theme_bw() +
  theme(plot.title = element_text(size = 9, vjust = 3, face = "italic"),
        plot.title.position = "plot",
        axis.title.y = element_text(size = 9),
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

# Fit quantile regression models: additive smooth term for year
mod_smoo_high <- gcrq(log_SDC ~ ps(log_fire_area, monotone = 1) + ps(Year),
                      lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_smoo_mix <- gcrq(log_SDC ~ ps(log_fire_area, monotone = 1) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_smoo_low <- gcrq(log_SDC ~ ps(log_fire_area, monotone = 1) + ps(Year),
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Fit quantile regression models: additive linear term for year
mod_lin_high <- gcrq(log_SDC ~ ps(log_fire_area, monotone = 1) + Year_std,
                     lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_lin_mix <- gcrq(log_SDC ~ ps(log_fire_area, monotone = 1) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_lin_low <- gcrq(log_SDC ~ ps(log_fire_area, monotone = 1) + Year_std,
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions for additive smooth
pred_smoo_high <- predict_year_fn(mod_smoo_high, fire_metrics_high)
pred_smoo_mix <- predict_year_fn(mod_smoo_mix, fire_metrics_mix)
pred_smoo_low <- predict_year_fn(mod_smoo_low, fire_metrics_low)

pred_smoo <- bind_rows(pred_smoo_high, pred_smoo_mix, pred_smoo_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Predictions for additive linear (with p-values)
pred_lin_high <- predict_year_fn(mod_lin_high, fire_metrics_high) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_high, "Year_std") ) # join p-values for linear terms
pred_lin_mix <- predict_year_fn(mod_lin_mix, fire_metrics_mix) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_mix, "Year_std") ) # join p-values for linear terms
pred_lin_low <- predict_year_fn(mod_lin_low, fire_metrics_low) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "tau", names_prefix = "fit_", values_to = "fit") %>%
  mutate(tau = as.numeric(tau)) %>%
  left_join( extract_p(mod_lin_low, "Year_std") ) # join p-values for linear terms

pred_lin <- bind_rows(pred_lin_high, pred_lin_mix, pred_lin_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))


# Plotting
p5 <- ggplot( filter(pred_smoo, log_fire_area == 3.5) ) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(limits = range(fire_metrics$log_SDC, na.rm = T), 
                     breaks = c(-3, -2.5, -2, -1.5),
                     labels = c("0.001", "0.003", "0.01", "0.03"),
                     name = "SDC\nparameter") +
  scale_x_continuous(name = "Year") +
  geom_ribbon(aes(x = Year, ymin = fit_0.05, ymax = fit_0.95, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = Year, ymin = fit_0.25, ymax = fit_0.75, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = Year, y = fit_0.05, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.25, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.75, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.95, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = Year, y = fit_0.5, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val >= 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.2, lty = 3) +
  geom_line(data = filter(pred_lin, log_fire_area == 3.5 & p_val < 0.05),
            mapping = aes(x = Year, y = fit, group = tau), lwd = 0.5, lty = 3) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")



# Combine plots ---------

plot_grid(p1, p2, p3, p4, p5, label_size = 10,
          label_y = c(0.85, 1, 1, 0.85, 1),
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)"),
          ncol = 1, align = "v", rel_heights = c(1.4, 1, 1, 1.2, 1.3))

# Export figure
ggsave(paste0("Figures/figure4_scaling_year.png"),
       width = 110, height = 150, units = "mm", dpi = 600)






