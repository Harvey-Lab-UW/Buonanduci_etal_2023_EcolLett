# Spatial scaling of burn severity

# Components of figure illustrating patch size and structure metrics

require(tidyverse)
require(sf)
require(raster)
require(R.utils)
require(cubature)
require(scales) # for numeric axis labels with commas


# Load function for fitting truncated lognormal distribution
source("Scripts/fit_truncated_ln.R")

# Function for plotting empirical patch distributions with
# fitted truncated lognormal distribution curve 
plot_truncated_ln <- function(patch_dist, xmin, xmax, beta, psi){
  
  # Filter patch size distribution to only include values >= xmin
  patch_dist <- patch_dist[patch_dist >= xmin]
  
  # Calculate normalization constant, a
  f_a <- function(x) { exp( -beta*log(x) - psi*log(x)^2 ) }
  a <- 1 / cubintegrate(f = f_a, lower = xmin, upper = xmax, method = "pcubature")$integral
  
  # Specify fitted lognormal density function
  f_ln <- function(x) { a * exp( -beta*log(x) - psi*log(x)^2 ) }
  
  # Loop through and calculate inverse cumulative probabilities
  # i.e., P(X >= x)
  x_range <- exp( seq(from = log(xmin), to = log(max(patch_dist)), by = 0.1) )
  P_predict <- c()
  
  for (i in 1:length(x_range)) {
    P <- cubintegrate(f = f_ln, lower = x_range[i], upper = xmax, method = "pcubature")$integral
    P_predict <- c(P_predict, P)
  }
  
  # Store predictions in tibble for plotting
  ln_predict <- tibble(x = x_range, P_x = P_predict)
  
  # Calculate inverse empirical cumulative density
  x_n <- length(patch_dist)
  x_ecdf <- c()
  for(x in 1:x_n){
    x_ecdf <- c(x_ecdf,
                sum(patch_dist >= patch_dist[x]) / x_n)
  }
  dist <- tibble(patch_area_ha = patch_dist, inv_ecdf = x_ecdf)
  
  # Plot inverse empirical cumulative density and fitted curve
  ggplot(data = dist, aes(x = patch_area_ha, y = inv_ecdf)) +
    geom_point(shape = 1, alpha = 0.6) +
    geom_line(data = ln_predict, aes(x = x, y = P_x), color = "#8FD744FF") +
    scale_x_log10(name = "Patch size (ha)") +
    scale_y_log10() 
}

# Function to calculate area-weighted means
AW_mean <- function(x){sum( x * (x / sum(x)))}



# Load data --------------

# Load the fire perimeters
fire_perims <- st_read("Data/fire_perims.gpkg") %>%
  mutate(fire_area = st_area(.)) %>%
  mutate(fire_area_ha = as.numeric(fire_area / 10000)) 

# Load fire metrics data 
fire_metrics <- read_csv("Data/fire_metrics.csv")

# Filter to demo fires for figure
fire_perims_demo <- bind_rows( filter(fire_perims, Fire_ID == "WY4437710988020190902"),
                             filter(fire_perims, Fire_ID == "OR4317212252420090913"),
                             filter(fire_perims, Fire_ID == "OR4310612253920020714") )


# Plot correlation of beta and psi parameters -------

ggplot(fire_metrics, aes(psi, beta)) +
  geom_point(color = "gray", alpha = 0.5) +
  geom_point(data = filter(fire_metrics, Fire_ID == "WY4437710988020190902"),
             mapping = aes(psi, beta), shape = 16, color = "#57A900") +
  geom_point(data = filter(fire_metrics, Fire_ID == "OR4317212252420090913"),
             mapping = aes(psi, beta), shape = 16, color = "#57A900") +
  geom_point(data = filter(fire_metrics, Fire_ID == "OR4310612253920020714"),
             mapping = aes(psi, beta), shape = 16, color = "#57A900") +
  labs(x = "\u03C8", y = "\u03B2") +
  theme_bw() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("Figures/figure2_beta_psi_correlation.png",
       width = 4, height = 1.5, units = "in", dpi = 600)




# Read in and plot high-severity rasters ----------

plot_raster_HS <- list()

for (i in 1:nrow(fire_perims_demo)){
  # Extract fire perimeter
  fire <- fire_perims_demo[i, ]
  
  # Read in raster
  high_sev <- raster( paste0("Data/tif/", fire$Fire_ID, "_HS.tif") ) %>%
    crop( extent(fire) ) %>% mask(fire)
  
  # Convert to data frame
  high_sev_df <- tibble(X_m = coordinates(high_sev)[,1],
                             Y_m = coordinates(high_sev)[,2],
                             severity = values(high_sev)) %>%
    drop_na(severity) %>%
    mutate(X_m = X_m - min(X_m)) %>%
    mutate(Y_m = Y_m - min(Y_m))
  
  # Plot
  plot_raster_HS[[i]] <- ggplot(high_sev_df, aes(x = X_m, y = Y_m, fill = factor(severity))) +
    geom_raster() +
    xlim(c(0, 15000)) + ylim(c(0, 15000)) +
    coord_fixed() +
    scale_fill_manual(
      values=c("#ECD5C0", "#873e23"),
      name = "Burn severity",
      labels = c("Low/moderate", "High")) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
}

# Save to files
plot_raster_HS[[1]]
ggsave(paste0("Figures/figure2_demo1_HS.png"),
       width = 2.7, height = 2.7, units = "in", dpi = 600)

plot_raster_HS[[2]]
ggsave(paste0("Figures/figure2_demo2_HS.png"),
       width = 2.7, height = 2.7, units = "in", dpi = 600)

plot_raster_HS[[3]] +
  theme(legend.position = c(0.75, 0.3),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.background=element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(paste0("Figures/figure2_demo3_HS.png"),
       width = 2.7, height = 2.7, units = "in", dpi = 600)




# Read in and plot distance-to-seed rasters ------------

plot_raster_DTS <- list()

for (i in 1:nrow(fire_perims_demo)){
  # Extract fire perimeter
  fire <- fire_perims_demo[i, ]
  
  # Read in raster
  distance_seed <- raster( paste0("Data/tif/", fire$Fire_ID, "_DTS.tif") ) %>%
    crop( extent(fire) ) %>% mask(fire)
  
  # Convert to data frame points
  distance_seed_df <- tibble(X_m = coordinates(distance_seed)[,1],
                             Y_m = coordinates(distance_seed)[,2],
                             DTS = values(distance_seed)) %>%
    drop_na(DTS) %>%
    mutate(X_m = X_m - min(X_m)) %>%
    mutate(Y_m = Y_m - min(Y_m))
  
  # Plot rasters
  plot_raster_DTS[[i]] <- ggplot(distance_seed_df, aes(x = X_m, y = Y_m, fill = DTS)) +
    geom_raster() +
    xlim(c(0, 15000)) + ylim(c(0, 15000)) +
    coord_fixed() +
    scale_fill_viridis_c(name = "Distance\nto seed (m)", option = "mako",
                         breaks = c(0, 500, 1000), limits = c(0, 1000), labels = comma) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
}


# Save to files
plot_raster_DTS[[1]]
ggsave(paste0("Figures/figure2_demo1_DTS.png"),
       width = 2.7, height = 2.7, units = "in", dpi = 600)

plot_raster_DTS[[2]]
ggsave(paste0("Figures/figure2_demo2_DTS.png"),
       width = 2.7, height = 2.7, units = "in", dpi = 600)

plot_raster_DTS[[3]] +
  theme(legend.position = c(0.75, 0.3),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.background=element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(paste0("Figures/figure2_demo3_DTS.png"),
       width = 2.7, height = 2.7, units = "in", dpi = 600)






# Patch size distribution plots ------------------

# Load the patch distributions
patch_dist <- read_csv("Data/fire_patch_distributions.csv")

# Loop through fire perimeters, fit truncated lognormal,
# extract beta and psi
ln_estimates <- list()
plot_patch_fits <- list()

for (i in 1:nrow(fire_perims_demo)){

  # Pull out fire
  fire <- fire_perims_demo[i, ]
  
  # Filter patch size distributions
  dist <- filter(patch_dist, Fire_ID == fire$Fire_ID) 
  
  # Check whether there are at least n=10 patches >=1 ha in size
  if (nrow(filter(dist, patch_area_ha >= 1)) >= 10) {
    
    # Check whether truncated lognormal MLE algorithm
    # (a) fails to converge within 5 minutes or (b) returns error
    t <- try( withTimeout(fit_truncated_ln(patch_dist = dist$patch_area_ha, 
                                           xmin = 1, xmax = fire$fire_area_ha, eta = 0.1),
                          timeout = 60*5), silent = TRUE)
    
    if ("try-error" %in% class(t)){
      ln_estimates[[i]] <- tibble(Fire_ID = fire$Fire_ID,
                                  Fire_Name = fire$Fire_Name,
                                  fire_area_ha = fire$fire_area_ha,
                                  beta = NA, psi = NA,
                                  n = nrow(dist),
                                  n_gte_xmin = sum(dist$patch_area_ha >= 1))
    } else {
      
      # If no error returned, fit truncated lognormal
      ln_mod <- fit_truncated_ln(patch_dist = dist$patch_area_ha, 
                                 xmin = 1, xmax = fire$fire_area_ha, eta = 0.1)
      
      # Add parameter estimates to list
      ln_estimates[[i]] <- tibble(Fire_ID = fire$Fire_ID,
                                  Fire_Name = fire$Fire_Name,
                                  fire_area_ha = fire$fire_area_ha,
                                  patch_AW_mean = AW_mean(dist$patch_area_ha),
                                  beta = ln_mod$beta,
                                  psi = ln_mod$psi,
                                  n = nrow(dist),
                                  n_gte_xmin = sum(dist$patch_area_ha >= 1))
      
      # Plot fit
      patch_vline <- tibble(x = c(AW_mean(dist$patch_area_ha), AW_mean(dist$patch_area_ha)),
                            y = c(8e-3, 0.4))
      
      plot_patch_fits[[i]] <- plot_truncated_ln(patch_dist = dist$patch_area_ha, 
                                                xmin = 1, xmax = fire$fire_area_ha,
                                                beta = ln_mod$beta, psi = ln_mod$psi) +
        geom_line(data = patch_vline, mapping = aes(x, y), lty = 3) +
        scale_x_log10(breaks = c(1, 100, 10000), labels = comma,
                      name = "Patch size (ha)") +
        coord_cartesian(xlim = c(1, 15000), ylim = c(7e-3, 1)) +
        theme_bw() +
        theme(axis.title = element_text(size = 8),
              axis.text = element_text(size = 8),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
    }
    
    
  } else {
    ln_estimates[[i]] <- tibble(Fire_ID = fire$Fire_ID,
                                Fire_Name = fire$Fire_Name,
                                fire_area_ha = fire$fire_area_ha,
                                beta = NA, psi = NA,
                                n = nrow(dist), 
                                n_gte_xmin = sum(dist$patch_area_ha >= 1))
  }
}

# Bind together list
ln_estimates <- bind_rows(ln_estimates)

plot_patch_fits[[1]] + ylab("Inverse cumulative\nprobability")
ggsave(paste0("Figures/figure2_demo1_patch_dist.png"),
       width = 2.1, height = 1.6, units = "in", dpi = 600)

plot_patch_fits[[2]] + ylab(" \n ")
ggsave(paste0("Figures/figure2_demo2_patch_dist.png"),
       width = 2.1, height = 1.6, units = "in", dpi = 600)

plot_patch_fits[[3]] + ylab(" \n ")
ggsave(paste0("Figures/figure2_demo3_patch_dist.png"),
       width = 2.1, height = 1.6, units = "in", dpi = 600)






# Distance to seed distribution plots ----------------

DTS_distributions <- list()
plot_DTS_fits <- list()

for (i in 1:nrow(fire_perims_demo)){
  # Pull out fire
  fire <- fire_perims_demo[i, ]
  
  # Distance to seed distribution
  DTS_dist <- read_csv("Data/fire_DTS_distributions.csv") %>%
    filter(Fire_ID == fire$Fire_ID) %>%
    mutate(DTS_area_ha = DTS_cells * 30 * 30 / 10000)
  
  # DTS distribution parameter
  DTS_dist_param <- fire_metrics %>%
    filter(Fire_ID == fire$Fire_ID)
  
  # DTS model predictions
  DTS_dist <- DTS_dist %>% 
    mutate(pred = 1 / (10^(DTS_dist_param$SDC * DTS_dist_m)) )
  
  # Add distribution to list
  DTS_distributions[[i]] <- DTS_dist
  
  # Dotted lines
  core <- DTS_dist$DTS_cells_prp[DTS_dist$DTS_dist_m==150]
  max_dist <- max(DTS_dist$DTS_dist_m)
  DTS_vline <- tibble(x = c(150, 150),
                      y = c(-1, core))
  DTS_hline <- tibble(x = c(150, 1000), # set to 1000 m instead of max_dist
                      y = rep(core, 2))
  
  plot_DTS_fits[[i]] <- ggplot(DTS_dist) +
    geom_point(mapping = aes(x = DTS_dist_m, y = DTS_cells_prp), shape = 1, alpha = 0.6) +
    geom_line(mapping = aes(x = DTS_dist_m, y = pred), color = "#8FD744FF") +
    geom_line(data = DTS_vline, mapping = aes(x, y), lty = 3) +
    geom_line(data = DTS_hline, mapping = aes(x, y), lty = 3) +
    scale_x_continuous(breaks = c(0, 500, 1000), labels = comma) +
    labs(x = "Distance to seed (m)", y = "Inverse cumulative\nproportion") +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 1000)) +
    theme_bw() +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}


plot_DTS_fits[[1]] +
  scale_y_continuous(sec.axis = sec_axis(~.*DTS_distributions[[1]]$DTS_area_ha[1],
                                         name = " ",
                                         labels = comma))
ggsave(paste0("Figures/figure2_demo1_DTS_dist.png"),
       width = 2.35, height = 1.6, units = "in", dpi = 600)

plot_DTS_fits[[2]] +
  scale_y_continuous(sec.axis = sec_axis(~.*DTS_distributions[[2]]$DTS_area_ha[1],
                                         name = " ",
                                         labels = comma)) + ylab(" ")
ggsave(paste0("Figures/figure2_demo2_DTS_dist.png"),
       width = 2.25, height = 1.6, units = "in", dpi = 600)

plot_DTS_fits[[3]] +
  scale_y_continuous(sec.axis = sec_axis(~.*DTS_distributions[[3]]$DTS_area_ha[1],
                                         name = "Inverse cumulative\narea (ha)",
                                         labels = comma)) + ylab(" ")
ggsave(paste0("Figures/figure2_demo3_DTS_dist.png"),
       width = 2.28, height = 1.6, units = "in", dpi = 600)






