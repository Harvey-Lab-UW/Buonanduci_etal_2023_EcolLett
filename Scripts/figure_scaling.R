# Scaling analysis
# Figure illustrating quantile regression scaling relationships

# Load packages
require(tidyverse)
require(quantreg)
require(quantregGrowth)
require(cowplot)

# Load data--------------

# Load landscape metrics data file
fire_metrics <- read_csv("Data/fire_metrics.csv")

# Create factors
fire_metrics <- fire_metrics %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

# Pull out unique fire regimes
fire_metrics_high <- filter(fire_metrics, Fire_Regime == "High")
fire_metrics_mix <- filter(fire_metrics, Fire_Regime == "Mixed")
fire_metrics_low <- filter(fire_metrics, Fire_Regime == "Low")

# Create facet labels
FRG_names <- c(
  `High` = "High severity regime",
  `Mixed` = "Mixed severity regime",
  `Low` = "Low severity regime"
)



# Area-wtd mean patch size ----

# Fit quantile regression models
mod_high <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Make predictions
pred_high <- fire_metrics_high %>% 
  drop_na(log_patch_area_AW_mean) %>% bind_cols( as_tibble(mod_high$fitted.values) )
pred_mix <- fire_metrics_mix %>% 
  drop_na(log_patch_area_AW_mean) %>% bind_cols( as_tibble(mod_mix$fitted.values) )
pred_low <- fire_metrics_low %>% 
  drop_na(log_patch_area_AW_mean) %>% bind_cols( as_tibble(mod_low$fitted.values) )

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))


# Plotting
p1a <- ggplot(pred) +
  facet_grid( ~ Fire_Regime, labeller = as_labeller(FRG_names)) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Area-weighted\nmean (ha)") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_point(aes(x = log_fire_area, y = log_patch_area_AW_mean, color = Fire_Regime), size = 0.6, alpha = 0.1) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.25`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.75`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  ggtitle("Patch size metrics") +
  theme_bw() +
  theme(plot.title = element_text(size = 9, vjust = 3, face = "italic"),
        plot.title.position = "plot",
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(9.5, 5.5, 5.5, 11.5), "points"))

p1b <- ggplot(pred) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Area-weighted\nmean (ha)") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  ggtitle(" ") +
  theme_bw() +
  theme(plot.title = element_text(size = 9, vjust = 3),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(9.5, 5.5, 5.5, 11.5), "points"))


# Patch distribution: beta ----

# Fit quantile regression models
mod_high <- gcrq(beta ~ ps(log_fire_area),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(beta ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(beta ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Make predictions
pred_high <- fire_metrics_high %>% drop_na(beta) %>% 
  bind_cols( as_tibble(mod_high$fitted.values) ) 
pred_mix <- fire_metrics_mix %>% drop_na(beta) %>% 
  bind_cols( as_tibble(mod_mix$fitted.values) ) 
pred_low <- fire_metrics_low %>% drop_na(beta) %>% 
  bind_cols( as_tibble(mod_low$fitted.values) ) 

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Plotting
p2a <- ggplot(pred) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(name = "\u03B2\nparameter") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_point(aes(x = log_fire_area, y = beta, color = Fire_Regime), size = 0.6, alpha = 0.1) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.25`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.75`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

p2b <- ggplot(pred) +
  scale_y_continuous(name = "\u03B2\nparameter") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# Patch distribution: psi ----

# Fit quantile regression models
mod_high <- gcrq(psi ~ ps(log_fire_area),
                    lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(psi ~ ps(log_fire_area),
                   lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(psi ~ ps(log_fire_area),
                   lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Make predictions
pred_high <- fire_metrics_high %>% drop_na(psi) %>% 
  bind_cols( as_tibble(mod_high$fitted.values) ) 
pred_mix <- fire_metrics_mix %>% drop_na(psi) %>% 
  bind_cols( as_tibble(mod_mix$fitted.values) ) 
pred_low <- fire_metrics_low %>% drop_na(psi) %>% 
  bind_cols( as_tibble(mod_low$fitted.values) ) 

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Plotting
p3a <- ggplot(pred) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(limits = c(-0.21, 1),
                     name = "\u03C8\nparameter") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_point(aes(x = log_fire_area, y = psi, color = Fire_Regime), size = 0.6, alpha = 0.1) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.25`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.75`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

p3b <- ggplot(pred) +
  scale_y_continuous(#limits = c(-0.21, 1),
    name = "\u03C8\nparameter") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# Core area ----

# Fit quantile regression models
mod_high <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Make predictions
pred_high <- fire_metrics_high %>% 
  drop_na(log_total_core) %>% bind_cols( as_tibble(mod_high$fitted.values) )
pred_mix <- fire_metrics_mix %>% 
  drop_na(log_total_core) %>% bind_cols( as_tibble(mod_mix$fitted.values) )
pred_low <- fire_metrics_low %>% 
  drop_na(log_total_core) %>% bind_cols( as_tibble(mod_low$fitted.values) )

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Plotting
p4a <- ggplot(pred) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Total core\narea (ha)") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_point(aes(x = log_fire_area, y = log_total_core, color = Fire_Regime), size = 0.6, alpha = 0.1) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.25`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.75`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  ggtitle("Patch structure metrics") +
  theme_bw() +
  theme(plot.title = element_text(size = 9, vjust = 3, face = "italic"),
        plot.title.position = "plot",
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

p4b <- ggplot(pred) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Total core\narea (ha)") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  ggtitle(" ") +
  theme_bw() +
  theme(plot.title = element_text(size = 9, vjust = 3),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




# Seed decay coefficient ----

# Fit quantile regression models
mod_high <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Make predictions
pred_high <- fire_metrics_high %>% 
  drop_na(log_SDC) %>% bind_cols( as_tibble(mod_high$fitted.values) )
pred_mix <- fire_metrics_mix %>% 
  drop_na(log_SDC) %>% bind_cols( as_tibble(mod_mix$fitted.values) )
pred_low <- fire_metrics_low %>% 
  drop_na(log_SDC) %>% bind_cols( as_tibble(mod_low$fitted.values) )

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Plotting
p5a <- ggplot(pred) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(breaks = c(-3, -2.5, -2, -1.5),
                     labels = c("0.001", "0.003", "0.01", "0.03"),
                     name = "SDC\nparameter") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  geom_point(aes(x = log_fire_area, y = log_SDC, color = Fire_Regime), size = 0.6, alpha = 0.1) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.25`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.75`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8,
                       name = "Historical fire regime", labels = c("Low", "Mixed", "High")) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8, 
                        name = "Historical fire regime", labels = c("Low", "Mixed", "High")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

p5b <- ggplot(pred) +
  scale_y_continuous(breaks = c(-3, -2.5, -2, -1.5),
                     labels = c("0.001", "0.003", "0.01", "0.03"),
                     name = "SDC\nparameter") +
  scale_x_continuous(limits = c(2.6, 5.64), breaks = c(3, 4, 5), 
                     labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`, fill = Fire_Regime), alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.5`, color = Fire_Regime)) +
  geom_line(aes(x = log_fire_area, y = `0.05`, color = Fire_Regime), lwd = 0.1) +
  geom_line(aes(x = log_fire_area, y = `0.95`, color = Fire_Regime), lwd = 0.1) +
  scale_fill_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  scale_color_viridis_d(direction = -1, option = "inferno", end = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




# Combine plots ---------
pgrid1 <- plot_grid(p1a, p2a, p3a, p4a, p5a, 
                   ncol = 1, align = "v", label_size = 10,
                   rel_heights = c( 1.3, 1, 1, 1.3, 1.3),
                   label_y = c(0.85, 1, 1, 0.85, 1),
                   labels = c("(a)", "(c)", "(e)", "(g)", "(i)"))

pgrid2 <- plot_grid(p1b, p2b, p3b,p4b, p5b,
                    ncol = 1, align = "v", label_size = 10,
                    rel_heights = c(1.3, 1, 1, 1.3, 1.3), 
                    label_y = c(0.85, 1, 1, 0.85, 1),
                    labels = c("(b)", "(d)", "(f)", "(h)", "(j)"))

pgrid <- plot_grid(pgrid1, pgrid2, ncol = 2, rel_widths = c(3,1.7))

# Extract legend that is laid out horizontally
legend_b <- get_legend(
  p5a + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
)

# Add the legend underneath
plot_grid(pgrid, legend_b, ncol = 1, rel_heights = c(1, 0.07),
          align = "v", axis = "l")

# Export figure
ggsave("Figures/figure3_scaling.png", bg = "white",
       width = 6.5, height = 7, units = "in", dpi = 600)



