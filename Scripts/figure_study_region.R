# Study area figure

require(tidyverse)
require(sf)
require(raster)


# Load data layers

# All fires
fire_perims <- st_read("Data/fire_perims.gpkg") %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 500) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("High", "Mixed", "Low")))

# US states
states <- st_read("Data/states.gpkg")

# Forested regions
region_NR <- st_read("Data/region_NR.gpkg") 
region_PNW <- st_read("Data/region_PNW.gpkg") 


# Plotting

# Plot extent
plt_extent <- extent(fire_perims)

ggplot() +
  xlim( c(plt_extent[1], plt_extent[2]) ) +
  ylim( c(plt_extent[3], plt_extent[4]) ) +
  labs( x = NULL, y = NULL ) +
  geom_sf(data = states, fill = "gray92", color = NA) +
  geom_sf(data = region_NR, fill = "tan", color = NA) +
  geom_sf(data = region_PNW, fill = "gray", color = NA) +
  geom_sf(data = states, fill = NA, color = "white") +
  geom_sf(data = fire_perims, mapping = aes(fill = Fire_Regime), color = NA) +
  scale_fill_viridis_d(option = "inferno",
                       name = "Historical fire regime",
                       labels = c("High-severity",
                                  "Mixed-severity",
                                  "Low-severity")) +
  annotate("text", 
           x = plt_extent[1] + 0.32*(plt_extent[2] - plt_extent[1]), 
           y = plt_extent[3] + 0.35*(plt_extent[4] - plt_extent[3]), 
           hjust = 0, color = "gray20", size = 3, fontface = "bold",
           label = "Pacific\nNorthwest") +
  annotate("text", 
           x = plt_extent[1] + 0.88*(plt_extent[2] - plt_extent[1]), 
           y = plt_extent[3] + 0.75*(plt_extent[4] - plt_extent[3]), 
           hjust = 0, color = "gray20", size = 3, fontface = "bold",
           label = "Northern\nRockies") +
  theme_bw() +
  theme(legend.position = c(0.6, 0.12),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(3, 'mm'),
        legend.background = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

ggsave(paste0("Figures/figure1_study_region.pdf"),
       width = 110, height = 95, units = "mm", dpi = 600)

