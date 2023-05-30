# Scaling analysis
# Evaluate whether relationship btwn fire size and spatial metrics
# varies by time period

# Load packages
require(tidyverse)
require(quantreg)
require(quantregGrowth)

# Load functions
# K-fold validation function
source("Scripts/kfold_fun.R")


# Load data--------------

# Load landscape metrics data file
fire_metrics <- read_csv("Data/fire_metrics.csv")

# Create factors
fire_metrics <- fire_metrics %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) %>%
  mutate(Time_Period = factor(Time_Period)) 

# Mean-center and standardize year for use as linear predictor
fire_metrics <- fire_metrics %>%
  mutate(Year_std = (Year - mean(Year))/sd(Year))

# Pull out unique fire regimes
fire_metrics_high <- filter(fire_metrics, Fire_Regime == "High")
fire_metrics_mix <- filter(fire_metrics, Fire_Regime == "Mixed")
fire_metrics_low <- filter(fire_metrics, Fire_Regime == "Low")

# Quantiles to evaluate
taus <- c(0.05,0.5,0.95)



# CROSS-VALIDATION -------------

# First, partition data into 10 subsamples
# Stratify subsamples so that each subsample contains the same proportion
#   of Early versus Late data points as in the original dataset

# Function to subsample data in a stratified manner
# Returns string of subsample indices for input data frame
subsample_fn <- function(df, n_s = 10) {
  # Get indices of each category
  ind_df_cat1 <- which(df$Time_Period == "Early")
  ind_df_cat2 <- which(df$Time_Period == "Late")
  
  # Get sample sizes for each category
  n_cat1 <- cut(c(1:length(ind_df_cat1)), n_s) %>%
    as.numeric() %>% table() %>% as.numeric()
  n_cat2 <- cut(c(1:length(ind_df_cat2)), n_s) %>%
    as.numeric() %>% table() %>% as.numeric()
  
  # Assign indices to subsample groups
  ind_df <- c()
  ind_ss <- c()
  
  # Loop through subsamples 1:9
  for (i in 1:(n_s-1)){
    # Randomly select dataframe indices
    samp_cat1 <- sample(ind_df_cat1, size = n_cat1[i], replace = FALSE)
    samp_cat2 <- sample(ind_df_cat2, size = n_cat2[i], replace = FALSE)
    # Add selected indices to list
    ind_df <- c(ind_df, samp_cat1, samp_cat2)
    ind_ss <- c(ind_ss, rep(i, length = n_cat1[i] + n_cat2[i]))
    # Remove selected indices from options for next iteration
    ind_df_cat1 <- ind_df_cat1[-which(ind_df_cat1 %in% samp_cat1)]
    ind_df_cat2 <- ind_df_cat2[-which(ind_df_cat2 %in% samp_cat2)]
  }
  
  # For subsample 10, add remaining indices
  ind_df <- c(ind_df, ind_df_cat1, ind_df_cat2)
  ind_ss <- c(ind_ss, rep(n_s, length = length(ind_df_cat1) + length(ind_df_cat2)))
  
  # Create tibble
  ind_tibble <- tibble(ind_df = ind_df, ind_ss = ind_ss) %>%
    arrange(ind_df)
  
  # Return subsample indices
  return(ind_tibble$ind_ss)
}


# For each FOLD:
# Fit candidate models to training dataset (dataset excluding 10% subsample)
# Calculate quantile loss function for testing dataset (10% subsample)

# Create list to summarize all metrics
kfold_val_sum <- list()



# Area-weighted mean patch size --------

# High severity
set.seed(12)

fire_metrics_high_ss <- fire_metrics_high %>%
  drop_na(log_patch_area_AW_mean) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_high_null <- kfold_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                           tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add <- kfold_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + ps(Year),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add_lin <- kfold_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + Year_std,
                              tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int <- kfold_fun(form = log_patch_area_AW_mean ~ Time_Period + 
                            ps(log_fire_area, monotone = 1, shared.pen = TRUE, by = Time_Period),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int_2d <- kfold_2d_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1), tau = taus, 
                                data1 = filter(fire_metrics_high_ss, Time_Period == "Early"),
                                data2 = filter(fire_metrics_high_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Mixed severity
set.seed(34)

fire_metrics_mix_ss <- fire_metrics_mix %>%
  drop_na(log_patch_area_AW_mean) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_mix_null <- kfold_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                          tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add <- kfold_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + ps(Year),
                         tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add_lin <- kfold_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + Year_std,
                             tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int <- kfold_fun(form = log_patch_area_AW_mean ~ Time_Period + 
                           ps(log_fire_area, monotone = 1, shared.pen = TRUE, by = Time_Period),
                         tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int_2d <- kfold_2d_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1), tau = taus, 
                               data1 = filter(fire_metrics_mix_ss, Time_Period == "Early"),
                               data2 = filter(fire_metrics_mix_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Low severity
set.seed(2345)

fire_metrics_low_ss <- fire_metrics_low %>%
  drop_na(log_patch_area_AW_mean) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_low_null <- kfold_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                          tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add <- kfold_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + ps(Year),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add_lin <- kfold_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1) + Year_std,
                             tau = taus, data = fire_metrics_low_ss) %>%
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add_lin") %>%
  group_by(resp, tau, model, Region, Fire_Regime) %>%
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int <- kfold_fun(form = log_patch_area_AW_mean ~ Time_Period + 
                           ps(log_fire_area, monotone = 1, shared.pen = TRUE, by = Time_Period),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int_2d <- kfold_2d_fun(form = log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1), tau = taus, 
                               data1 = filter(fire_metrics_low_ss, Time_Period == "Early"),
                               data2 = filter(fire_metrics_low_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))



# Add results to list
kfold_val_sum <- bind_rows(kfold_val_sum, 
                           val_high_null, val_high_add, val_high_int, val_high_int_2d,
                           val_mix_null, val_mix_add, val_mix_int, val_mix_int_2d,
                           val_low_null, val_low_add, val_low_int, val_low_int_2d)



# Patch size distribution: beta --------

# High severity
set.seed(129)

fire_metrics_high_ss <- fire_metrics_high %>%
  drop_na(beta) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_high_null <- kfold_fun(form = beta ~ ps(log_fire_area),
                           tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add <- kfold_fun(form = beta ~ ps(log_fire_area) + ps(Year),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add_lin <- kfold_fun(form = beta ~ ps(log_fire_area) + Year_std,
                              tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int <- kfold_fun(form = beta ~ Time_Period + 
                            ps(log_fire_area, shared.pen = TRUE, by = Time_Period),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int_2d <- kfold_2d_fun(form = beta ~ ps(log_fire_area), tau = taus, 
                                data1 = filter(fire_metrics_high_ss, Time_Period == "Early"),
                                data2 = filter(fire_metrics_high_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Mixed severity
set.seed(951)

fire_metrics_mix_ss <- fire_metrics_mix %>%
  drop_na(beta) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_mix_null <- kfold_fun(form = beta ~ ps(log_fire_area),
                          tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add <- kfold_fun(form = beta ~ ps(log_fire_area) + ps(Year),
                         tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add_lin <- kfold_fun(form = beta ~ ps(log_fire_area) + Year_std,
                             tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int <- kfold_fun(form = beta ~ Time_Period + 
                           ps(log_fire_area, shared.pen = TRUE, by = Time_Period),
                         tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int_2d <- kfold_2d_fun(form = beta ~ ps(log_fire_area), tau = taus, 
                               data1 = filter(fire_metrics_mix_ss, Time_Period == "Early"),
                               data2 = filter(fire_metrics_mix_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Low severity
set.seed(753)

fire_metrics_low_ss <- fire_metrics_low %>%
  drop_na(beta) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_low_null <- kfold_fun(form = beta ~ ps(log_fire_area),
                          tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add <- kfold_fun(form = beta ~ ps(log_fire_area) + ps(Year),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add_lin <- kfold_fun(form = beta ~ ps(log_fire_area) + Year_std,
                             tau = taus, data = fire_metrics_low_ss) %>%
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add_lin") %>%
  group_by(resp, tau, model, Region, Fire_Regime) %>%
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int <- kfold_fun(form = beta ~ Time_Period + 
                           ps(log_fire_area, shared.pen = TRUE, by = Time_Period),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int_2d <- kfold_2d_fun(form = beta ~ ps(log_fire_area), tau = taus, 
                               data1 = filter(fire_metrics_low_ss, Time_Period == "Early"),
                               data2 = filter(fire_metrics_low_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Add results to list
kfold_val_sum <- bind_rows(kfold_val_sum, 
                           val_high_null, val_high_add, val_high_int, val_high_int_2d,
                           val_mix_null, val_mix_add, val_mix_int, val_mix_int_2d,
                           val_low_null, val_low_add, val_low_int, val_low_int_2d)




# Patch size distribution: psi --------

# High severity
set.seed(129)

fire_metrics_high_ss <- fire_metrics_high %>%
  drop_na(psi) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_high_null <- kfold_fun(form = psi ~ ps(log_fire_area),
                           tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add <- kfold_fun(form = psi ~ ps(log_fire_area) + ps(Year),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add_lin <- kfold_fun(form = psi ~ ps(log_fire_area) + Year_std,
                              tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int <- kfold_fun(form = psi ~ Time_Period + 
                            ps(log_fire_area, shared.pen = TRUE, by = Time_Period),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int_2d <- kfold_2d_fun(form = psi ~ ps(log_fire_area), tau = taus, 
                                data1 = filter(fire_metrics_high_ss, Time_Period == "Early"),
                                data2 = filter(fire_metrics_high_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Mixed severity
set.seed(951)

fire_metrics_mix_ss <- fire_metrics_mix %>%
  drop_na(psi) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_mix_null <- kfold_fun(form = psi ~ ps(log_fire_area),
                          tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add <- kfold_fun(form = psi ~ ps(log_fire_area) + ps(Year),
                         tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add_lin <- kfold_fun(form = psi ~ ps(log_fire_area) + Year_std,
                             tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int <- kfold_fun(form = psi ~ Time_Period + 
                           ps(log_fire_area, shared.pen = TRUE, by = Time_Period),
                         tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int_2d <- kfold_2d_fun(form = psi ~ ps(log_fire_area), tau = taus, 
                               data1 = filter(fire_metrics_mix_ss, Time_Period == "Early"),
                               data2 = filter(fire_metrics_mix_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Low severity
set.seed(753)

fire_metrics_low_ss <- fire_metrics_low %>%
  drop_na(psi) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_low_null <- kfold_fun(form = psi ~ ps(log_fire_area),
                          tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add <- kfold_fun(form = psi ~ ps(log_fire_area) + ps(Year),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add_lin <- kfold_fun(form = psi ~ ps(log_fire_area) + Year_std,
                             tau = taus, data = fire_metrics_low_ss) %>%
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add_lin") %>%
  group_by(resp, tau, model, Region, Fire_Regime) %>%
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int <- kfold_fun(form = psi ~ Time_Period + 
                           ps(log_fire_area, shared.pen = TRUE, by = Time_Period),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int_2d <- kfold_2d_fun(form = psi ~ ps(log_fire_area), tau = taus, 
                               data1 = filter(fire_metrics_low_ss, Time_Period == "Early"),
                               data2 = filter(fire_metrics_low_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Add results to list
kfold_val_sum <- bind_rows(kfold_val_sum, 
                           val_high_null, val_high_add, val_high_int, val_high_int_2d,
                           val_mix_null, val_mix_add, val_mix_int, val_mix_int_2d,
                           val_low_null, val_low_add, val_low_int, val_low_int_2d)





# Core area --------

# High severity
set.seed(131)

fire_metrics_high_ss <- fire_metrics_high %>%
  drop_na(log_total_core) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_high_null <- kfold_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1),
                           tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add <- kfold_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1) + ps(Year),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add_lin <- kfold_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1) + Year_std,
                              tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int <- kfold_fun(form = log_total_core ~ Time_Period + 
                            ps(log_fire_area, monotone = 1, shared.pen = TRUE, by = Time_Period),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int_2d <- kfold_2d_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1), tau = taus, 
                                data1 = filter(fire_metrics_high_ss, Time_Period == "Early"),
                                data2 = filter(fire_metrics_high_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Mixed severity
set.seed(22)

fire_metrics_mix_ss <- fire_metrics_mix %>%
  drop_na(log_total_core) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_mix_null <- kfold_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1),
                           tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add <- kfold_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1) + ps(Year),
                          tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add_lin <- kfold_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1) + Year_std,
                              tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int <- kfold_fun(form = log_total_core ~ Time_Period + 
                            ps(log_fire_area, monotone = 1, shared.pen = TRUE, by = Time_Period),
                          tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int_2d <- kfold_2d_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1), tau = taus, 
                                data1 = filter(fire_metrics_mix_ss, Time_Period == "Early"),
                                data2 = filter(fire_metrics_mix_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Low severity
set.seed(129)

fire_metrics_low_ss <- fire_metrics_low %>%
  drop_na(log_total_core) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_low_null <- kfold_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1),
                          tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add <- kfold_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1) + ps(Year),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add_lin <- kfold_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1) + Year_std,
                             tau = taus, data = fire_metrics_low_ss) %>%
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add_lin") %>%
  group_by(resp, tau, model, Region, Fire_Regime) %>%
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int <- kfold_fun(form = log_total_core ~ Time_Period + 
                           ps(log_fire_area, monotone = 1, shared.pen = TRUE, by = Time_Period),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int_2d <- kfold_2d_fun(form = log_total_core ~ ps(log_fire_area, monotone = 1), tau = taus, 
                               data1 = filter(fire_metrics_low_ss, Time_Period == "Early"),
                               data2 = filter(fire_metrics_low_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Add results to list
kfold_val_sum <- bind_rows(kfold_val_sum, 
                           val_high_null, val_high_add, val_high_int, val_high_int_2d,
                           val_mix_null, val_mix_add, val_mix_int, val_mix_int_2d,
                           val_low_null, val_low_add, val_low_int, val_low_int_2d)




# Seed decay coefficient --------

# High severity
set.seed(129)

fire_metrics_high_ss <- fire_metrics_high %>%
  drop_na(log_SDC) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_high_null <- kfold_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1),
                           tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add <- kfold_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1) + ps(Year),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_add_lin <- kfold_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1) + Year_std,
                              tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int <- kfold_fun(form = log_SDC ~ Time_Period + 
                            ps(log_fire_area, monotone = -1, shared.pen = TRUE, by = Time_Period),
                          tau = taus, data = fire_metrics_high_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_high_int_2d <- kfold_2d_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1), tau = taus, 
                                data1 = filter(fire_metrics_high_ss, Time_Period == "Early"),
                                data2 = filter(fire_metrics_high_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "High", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Mixed severity
set.seed(951)

fire_metrics_mix_ss <- fire_metrics_mix %>%
  drop_na(log_SDC) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_mix_null <- kfold_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1),
                          tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add <- kfold_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1) + ps(Year),
                         tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_add_lin <- kfold_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1) + Year_std,
                             tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "add_lin") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int <- kfold_fun(form = log_SDC ~ Time_Period + 
                           ps(log_fire_area, monotone = -1, shared.pen = TRUE, by = Time_Period),
                         tau = taus, data = fire_metrics_mix_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_mix_int_2d <- kfold_2d_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1), tau = taus, 
                               data1 = filter(fire_metrics_mix_ss, Time_Period == "Early"),
                               data2 = filter(fire_metrics_mix_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Mixed", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Low severity
set.seed(753)

fire_metrics_low_ss <- fire_metrics_low %>%
  drop_na(log_SDC) %>% # drop missing response var
  mutate(subsample = subsample_fn(.))

val_low_null <- kfold_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1),
                          tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "null") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add <- kfold_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1) + ps(Year),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_add_lin <- kfold_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1) + Year_std,
                             tau = taus, data = fire_metrics_low_ss) %>%
  mutate(Region = "Combined", Fire_Regime = "Low", model = "add_lin") %>%
  group_by(resp, tau, model, Region, Fire_Regime) %>%
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int <- kfold_fun(form = log_SDC ~ Time_Period + 
                           ps(log_fire_area, monotone = -1, shared.pen = TRUE, by = Time_Period),
                         tau = taus, data = fire_metrics_low_ss) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))

val_low_int_2d <- kfold_2d_fun(form = log_SDC ~ ps(log_fire_area, monotone = -1), tau = taus, 
                               data1 = filter(fire_metrics_low_ss, Time_Period == "Early"),
                               data2 = filter(fire_metrics_low_ss, Time_Period == "Late")) %>% 
  mutate(Region = "Combined", Fire_Regime = "Low", model = "int_2d") %>% 
  group_by(resp, tau, model, Region, Fire_Regime) %>% 
  summarize(total_loss = sum(loss), avg_loss = mean(loss))


# Add results to list
kfold_val_sum <- bind_rows(kfold_val_sum, 
                           val_high_null, val_high_add, val_high_int, val_high_int_2d,
                           val_mix_null, val_mix_add, val_mix_int, val_mix_int_2d,
                           val_low_null, val_low_add, val_low_int, val_low_int_2d)



# Create summary table
kfold_val_table <- kfold_val_sum %>%
  group_by(resp, model, Fire_Regime) %>%
  summarise(avg_loss = mean(avg_loss))

# Write summary table to file
write_csv(kfold_val_table, "Data/cross_validation_time.csv")

