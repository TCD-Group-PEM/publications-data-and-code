
# create color blind friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(ggforce)
library(ggtext)
library(tidyverse)


############################################################################################
######################################## light part ########################################
############################################################################################


lightpart <- function(rate_of_change, current_daylength) {
  ifelse(current_daylength > 15, 0, pmax(pmin(rate_of_change / current_daylength, 1), 0))
}

calculate_day_length <- function(latitude, day_of_year) {
  # Convert latitude to radians
  latitude_rad <- latitude * pi / 180
  
  # Calculate declination angle for the given day of the year
  declination <- 23.45 * sin(2 * pi / 365 * (day_of_year - 81))
  
  # Calculate acos_arg
  acos_arg <- -tan(latitude_rad) * tan(declination * atan(1.0) * 4.0 / 180.0)
  
  # Initialize vector to store day length
  current_daylength <- numeric(length(acos_arg))
  
  # Check if argument is within the valid range for ACOS function
  current_daylength[acos_arg > 1.0] <- 0.0
  current_daylength[acos_arg < -1.0] <- 24.0
  
  # Calculate day length using the formula for the remaining values
  valid_indices <- acos_arg >= -1.0 & acos_arg <= 1.0
  current_daylength[valid_indices] <- (24.0 / atan(1.0) / 4.0) * (acos(acos_arg[valid_indices]) +
                                                                    (atan(1.0) * 0.0) / 180.0 * sin(latitude_rad) * sin(declination[valid_indices] * atan(1.0) * 4.0 / 180.0))
  
  return(current_daylength)
}

adjust_day_length <- function(daylength) {
  # Replace NA values in daylength vector with zeros
  daylength[is.na(daylength)] <- 0
  
  # Initialize a vector to store adjusted day lengths
  adjusted_day_lengths <- numeric(length(daylength))
  
  # Loop through the vector
  for (i in 2:length(daylength)) {
    # Subtract the current day length from the previous day length timestep value
    adjusted_day_lengths[i] <- (max(daylength[i-1] - daylength[i], 0))
  }
  
  return(adjusted_day_lengths)
}

equator      <- 0
low_latitude <- 23
mid_latitude  <- 45
high_latitude <- 67

latveclight  <- c(equator, low_latitude, mid_latitude, high_latitude)

day_of_year  <- 1:365

# calculate_day_length(latitude = , day_of_year = )
# adjust_day_length(daylength = )

dflight <- data.frame(
  doy =  rep(NA, 365),
  eq =   rep(NA, 365),
  low =  rep(NA, 365),
  mid =  rep(NA, 365),
  high = rep(NA, 365)
)

dflight[1] <- day_of_year
dflight[2] <- calculate_day_length(latitude = equator, day_of_year = day_of_year)
dflight[3] <- calculate_day_length(latitude = low_latitude, day_of_year = day_of_year)
dflight[4] <- calculate_day_length(latitude = mid_latitude, day_of_year = day_of_year)
dflight[5] <- calculate_day_length(latitude = high_latitude, day_of_year = day_of_year)


# List to store rate of change values
lightplot_roc <- list()

# Column names corresponding to latitudes in dflight
column_names <- c("eq", "low", "mid", "high")  

# Calculate adjusted day lengths for each latitude
for (i in seq_along(column_names)) {
  lightplot_roc[[i]] <- adjust_day_length(dflight[[column_names[i]]])
}

lightplot_roc_df <- as.data.frame(do.call(cbind, lightplot_roc))
colnames(lightplot_roc_df) <- c("eq_roc", "low_roc", "mid_roc", "high_roc")  # Rename columns

# Merge with dflight
dflight_merged <- cbind(dflight, lightplot_roc_df)

# Apply lightpart() to each corresponding pair of columns
dflight_roc <- dflight_merged %>%
  mutate(across(
    .cols = all_of(paste0(column_names, "_roc")),  # Select *_roc columns
    .fns = ~ lightpart(.x, dflight_merged[[sub("_roc$", "", cur_column())]]),  # Apply lightpart
    .names = "{.col}_lightpart"  # Rename new columns as *_roc_adjusted
  ))


############################################################################################
######################################### temp part ########################################
############################################################################################

tempparteq <- function(T_Air) {
  ifelse(T_Air < 0, 1, exp(- (0.5 * T_Air) ^ 2))
}


dftemp <- data.frame(
  T_Air = c(seq(from = 10 , to = -2, length.out = 30))
)

dftemp$temp_response <- tempparteq(T_Air = dftemp$T_Air)

############################################################################################
##################################### moisture part ########################################
############################################################################################

# moistpart(:) = 1.0_wp - MIN(1.0_wp, w_soil_root_pot(:) / lctlib_phi_leaf_min) 

moistparteq <- function(w_soil_root_pot){
  w_soil_root_pot / -1.55
}

dfmoist <- data.frame(
  w_soil_root_pot = c(seq(from = 0 , to = -1.55, length.out = 30))
)

dfmoist$moist_response <- moistparteq(w_soil_root_pot = dfmoist$w_soil_root_pot)


############################################################################################
################################## draw the Figure 1 #######################################
############################################################################################



lightplotsplit <- ggplot(data = dflight_roc)+
  geom_line( aes(x= doy, y = (eq_roc_lightpart), color = "eq"), lwd = 2) +
  geom_line( aes(x= doy, y = (low_roc_lightpart), color = "low"), lwd = 2) +
  geom_line( aes(x= doy, y = (mid_roc_lightpart), color = "mid"), lwd = 2) +
  geom_line( aes(x= doy, y = (high_roc_lightpart), color = "high"), lwd = 2) +
  scale_color_manual(
    values = c("eq" = cbbPalette[1], "low" = cbbPalette[2], "mid" = cbbPalette[3], "high" = cbbPalette[4]),  # Define color for the label
    labels = c("eq" = expression(0 * degree), "low" = expression(23 * degree), "mid" = expression(45 * degree), "high" = expression(67 * degree)),
    breaks = c("eq", "low", "mid", "high")
  ) +
  guides(
    color = guide_legend(order = 1)  # Make sure legend items appear in the specified order
  ) +
  labs(
    title = "Light component of the model at different latitudes throughout the year",
    x = "Day of Year",
    color = "Light response at °N lat.") +
  ylab(expression(f[turn]^{leaf}))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.x = element_markdown(size = 36, face = 'plain'),
    axis.title.y = element_text(size = 36),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 28, face = 'plain'),
    legend.title = element_text(size = 28, face = 'bold'),
    legend.position = "inside",
    legend.position.inside = c(0, .9),   # (x,y) in npc coordinates
    legend.justification = c("left", "top"),  # how the box anchors
    legend.background = element_rect(fill = "white", color = "black")  # optional
  )  +
  facet_zoom(ylim = c(-0.001, 0.02)); print(lightplotsplit)


tempplotsplit <- ggplot(data = dftemp)+
  geom_line( aes(x= T_Air, y = temp_response, color = "Temperature response"), lwd = 2) +
  labs(
    title = "Temperature component of the model",
    x = paste0( "Air temperature", "\t \u00B0", "C"),
    y = "") +
  scale_x_reverse() +
  ylab(expression(f[turn]^{leaf}))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.x = element_markdown(size = 36, face = 'plain'),
    axis.title.y = element_text(size = 36),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 28, face = 'plain'),
    legend.title = element_text(size = 28, face = 'bold'),
    legend.position = "inside",
    legend.position.inside = c(0, .9),   # (x,y) in npc coordinates
    legend.justification = c("left", "top"),  # how the box anchors
    legend.background = element_rect(fill = "white", color = "black")  # optional
  ) +
  scale_color_manual(
    name = "Data",
    breaks = c("Temperature response"),
    values = c("Temperature response" = 'red')
  ) ; print(tempplotsplit)


moistplotsplit <- ggplot(data = dfmoist)+
  geom_line( aes(x= w_soil_root_pot, y = moist_response, color = "Moisture response"), lwd = 2) +
  labs(
    title = "Moisture component of the model",
    x = "Soil Water Potential (MPa)",
    y = "") +
  ylab(expression(f[turn]^{leaf}))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.title.y = element_text(size = 36),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 28, face = 'plain'),
    legend.title = element_text(size = 28, face = 'bold'),
    legend.position = "inside",
    legend.position.inside = c(1, 1),   # (x,y) in npc coordinates
    legend.justification = c("right", "top"),  # how the box anchors
    legend.background = element_rect(fill = "white", color = "black")  # optional
  ) +
  scale_color_manual(
    name = "Data",
    breaks = c("Moisture response"),
    values = c("Moisture response" = 'blue')
  ) ; print(moistplotsplit)


allresponse_functions <- ggarrange (
  ggarrange(tempplotsplit, 
            moistplotsplit, 
            ncol = 2, 
            labels = c("a)", "b)"), 
            label.x = 0.195,
            label.y = 0.95,
            font.label = list(size = 36)),
  
  lightplotsplit, 
  nrow = 2, 
  labels = c("", "c)"),
  label.y = c(1, .95),  # "1" puts it at the top of each subpanel
  label.x = c(0.03, 0.11),
  font.label = list(size = 36)
) ; print(allresponse_functions)


ggsave(filename = "/path/to/figure.png", 
       plot = allresponse_functions, 
       scale = 6,
       width = 1000,
       height = 1000,
       units = "px")







