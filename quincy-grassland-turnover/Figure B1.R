

# load libraries
library(tidyverse)
library(maps)
library(ggthemes)



########################################################################
################################ SETUP #################################
########################################################################


# this uses the data frame (US_ICs_EC_daily) with US-Tol (US-ICs) EC flux GPP from file 'Figure 3 and Table 4.R: 
ustol_annual_sum_map_data <- US_ICs_EC_daily %>%
  rename(GPP =GPP_NT) %>%
  mutate(year = year(Date),
         Time = yday(Date),
         GPP = GPP * 12.0107) %>%
  group_by(year) %>%
  filter(n() >= 365) %>%
  mutate(GPP = ifelse(GPP < 0, 0, GPP)) %>%
  summarize(total_GPP = sum(GPP, na.rm = TRUE), .groups = "drop") %>%
  summarize(
    mean_annual_GPP = mean(total_GPP, na.rm = TRUE),
    sd_annual_GPP = sd(total_GPP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site = 'US_Tol',
         lat = 68.75,
         lon = -149.25 )


# this uses the data frame (maj_ec_flux_total) with ES-LMa EC flux GPP from file 'Figure 3 and Table 4.R: 
eslma_annual_sum_map_data <- maj_ec_flux_total %>%
  select(date, GPP) %>%
  mutate(year = year(date),
         Time = yday(date),
         GPP = GPP * 12.0107) %>%
  group_by(year) %>%
  filter(n() >= 365) %>%
  mutate(GPP = ifelse(GPP < 0, 0, GPP)) %>%
  summarize(total_GPP = sum(GPP, na.rm = TRUE), .groups = "drop") %>%
  summarize(
    mean_annual_GPP = mean(total_GPP, na.rm = TRUE),
    sd_annual_GPP = sd(total_GPP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site = 'ES_LMa',
         lat = 39.94,
         lon = -5.77 )


# this uses the data frames (combined_data, plumber2sites_edit) with all PLUMBER2 sites EC flux GPP and site info from file 'Figure 5, Figure B3 and Table 5, Table C2.R':
map_for_annual_gpp_data <- combined_data %>%
  mutate(year = year(date),
         GPP = GPP * 12.0107) %>%
  group_by(site, year) %>%
  filter(n() >= 365) %>%
  mutate(GPP = ifelse(GPP < 0, 0, GPP)) %>%
  summarize(total_GPP = sum(GPP, na.rm = TRUE), .groups = "drop") %>%
  group_by(site) %>%
  summarize(
    mean_annual_GPP = mean(total_GPP, na.rm = TRUE),
    sd_annual_GPP = sd(total_GPP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(plumber2sites_edit %>%
              select(Site.ID, lat, lon,koppen_class_name) %>%
              mutate(Site.ID = gsub("-", "_", Site.ID)) %>%
              rename(site = Site.ID),
            by = 'site') %>%
  bind_rows(eslma_annual_sum_map_data) %>%
  bind_rows(ustol_annual_sum_map_data)




# Load world map
world_map <- map_data("world")

# Create the map
annual_gpp_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "gray90", color = "white") +
  geom_point(data = map_for_annual_gpp_data, 
             aes(x = lon, y = lat, color = mean_annual_GPP), 
             size = 4, alpha = 0.8) + 
  theme_minimal() +  
  labs(title = "Locations of sites ", x = "Longitude", y = "Latitude") +
  theme(legend.position = "right") +
  scale_color_viridis_c(name = expression("Annual GPP (gC " * m^-2 * " year"^-1 * ")"), option = "plasma") +
  coord_cartesian(xlim = c(-150, 180), ylim = c(-50, 80))  ; print(annual_gpp_map)
