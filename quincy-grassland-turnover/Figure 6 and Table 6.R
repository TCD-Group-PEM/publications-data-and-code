

# libraries
library(tidyverse)
library(data.table)
library(ggtext)

########################################################################
################################ SETUP #################################
########################################################################

find_earliest_min_gpp_all <- function(df, senescence_dates, eos_date) {
  df <- df %>%
    mutate(year = year(date)) %>%  
    left_join(senescence_dates, by = "year") %>%  
    filter(date >= as.Date(eos_date)) %>%  
    group_by(year) %>%
    slice_min(order_by = rollmean(GPP, k = 14, fill = NA, align = "center"), with_ties = FALSE) %>%
    ungroup()
  
  return(df)
}

# load all data sources from file Figure 5 'Figure 5, Figure B3 and Table 5, Table C2.R':
# - 'matched_data' 
# - 'matched_data_dev' 
# - 'combined_data' 
# - 'all_quincy_results' 
# - 'all_quincy_results_mod' 
# - 'summary_df' 
# - 'plumber2sites_edit' 
# - 'filtered_long_data_frame' 


################################################################
################################################################


################################################################

# calculate mean senescence periods dynamic model

################################################################



filtered_long_data_frame2 <- filtered_long_data_frame %>%
  pivot_longer(
    cols = c(DER.sos, DER.eos, TRS2.sos, TRS2.eos),
    names_to = "metric",
    values_to = "value" 
  ) %>%
  select(site, origin, metric, value)  

summary_df_yearly <- filtered_long_data_frame2 %>%
  mutate(value = ifelse(value < 0, 365 + value, value),
         value = ifelse(value > 365, value - 365, value))


grass_new_eos_cleaned <- lapply(all_quincy_results_mod, function(site_list) {
  list(doy = site_list$doy$Elmore) 
})

grass_new_eos_cleaned <- map(grass_new_eos_cleaned, ~ .x$doy)

# drop unused columns
grass_new_eos_cleaned <- map(grass_new_eos_cleaned, ~ .x[, .(origin, DER.sos, DER.eos)])

# convert to long df
grass_new_eos_long <- rbindlist(
  lapply(names(grass_new_eos_cleaned), function(site_name) {
    site_data <- grass_new_eos_cleaned[[site_name]]
    site_data[, site := site_name] 
    melt(site_data, id.vars = c("site", "origin"), variable.name = "metric", value.name = "value")
  })
)

# summarize by site and metric
grass_new_eos_long_summary <- grass_new_eos_long %>%
  mutate(value = ifelse(value < 0, 365 + value, value),
         value = ifelse(value > 365, value - 365, value)) %>%
  group_by(site, metric) %>% 
  summarise(
    mean_eos = mean(value, na.rm = TRUE),
    sd_eos = sd(value, na.rm = TRUE),
    .groups = "drop" 
  )

grass_new_eos_long_summary <- grass_new_eos_long_summary %>%
  left_join(plumber2sites_edit[,c(1,4, 30)], by = c("site" = "Site.ID"))

grass_new_eos_long_summary_fin <- grass_new_eos_long_summary %>%
  inner_join(summary_df, by = c("site", "metric", "PFT", "koppen_class_name"), 
             suffix = c(".qm", ".ec"))



################################################################

# repeat for default mode

################################################################


grass_dev_eos_cleaned <- lapply(all_quincy_results, function(site_list) {
  list(doy = site_list$doy$Elmore) 
})

grass_dev_eos_cleaned <- map(grass_dev_eos_cleaned, ~ .x$doy)

grass_dev_eos_cleaned <- map(grass_dev_eos_cleaned, ~ .x[, .(origin, DER.sos, DER.eos)])

grass_dev_eos_long <- rbindlist(
  lapply(names(grass_dev_eos_cleaned), function(site_name) {
    site_data <- grass_dev_eos_cleaned[[site_name]]
    site_data[, site := site_name] 
    melt(site_data, id.vars = c("site", "origin"), variable.name = "metric", value.name = "value")
  })
)

grass_dev_eos_long_summary <- grass_dev_eos_long %>%
  mutate(value = ifelse(value < 0, 365 + value, value),
         value = ifelse(value > 365, value - 365, value)) %>%
  group_by(site, metric) %>%
  summarise(
    mean_eos = mean(value, na.rm = TRUE),
    sd_eos = sd(value, na.rm = TRUE),
    .groups = "drop" 
  )

grass_dev_eos_long_summary <- grass_dev_eos_long_summary %>%
  left_join(plumber2sites_edit[,c(1,4, 30)], by = c("site" = "Site.ID"))

grass_dev_eos_long_summary_fin <- grass_dev_eos_long_summary %>%
  inner_join(summary_df, by = c("site", "metric", "PFT", "koppen_class_name"), 
             suffix = c(".dev", ".ec"))


################################################################
# clean up dual growing seasons or no growing seasons
################################################################

# remove all dual growing season years for cleaner results
invalid_years <- grass_new_eos_long %>%
  filter(metric %in% c("DER.eos", "DER.sos")) %>%
  group_by(site, origin) %>%
  filter(!is.na(value)) %>%
  summarize(n = n(), .groups = "drop") %>%
  filter(n > 2) # get rid of multiple growing seasons

invalid_years_flux <- summary_df_yearly %>%
  filter(metric %in% c("DER.eos", "DER.sos")) %>%
  group_by(site, origin) %>%
  filter(!is.na(value)) %>%
  summarize(n = n(), .groups = "drop") %>%
  filter(n > 2)  

invalid_years_dev <- grass_dev_eos_long %>%
  filter(metric %in% c("DER.eos", "DER.sos")) %>%
  group_by(site, origin) %>%
  filter(!is.na(value)) %>%
  summarize(n = n(), .groups = "drop") %>%
  filter(n > 2)  

all_invalid_years <- bind_rows(
  invalid_years %>% mutate(source = "grass_new_eos_long"),
  invalid_years_dev %>% mutate(source = "grass_dev_eos_long")
) %>%
  distinct(site, origin)  

unique(all_invalid_years$site)

filtered_invalid_years <- all_invalid_years %>%
  anti_join(invalid_years_flux, by = c("site", "origin"))



####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

################################################################
# prep and extract senescence period from flux data
################################################################


flux_eos_mean_first_then_min <- summary_df_yearly %>% 
  filter(metric %in% c("DER.eos", "DER.sos")) %>%
  anti_join(invalid_years_flux, by = c("site", "origin")) %>%
  filter(!is.na(value)) %>%
  mutate(value = ifelse(value < 0, 365 + value, value), 
         value = ifelse(value > 365, value - 365, value)) %>% 
  pivot_wider(names_from = metric, values_from = value) %>%
  group_by(site) %>%
  summarize(mean_eos = mean(DER.eos, na.rm =T),
            sd_eos = sd(DER.eos, na.rm =T),
            mean_sos = mean(DER.sos, na.rm =T),
            sd_sos = sd(DER.sos, na.rm =T)) %>%
  ungroup()


flux_eos_gpp_mean_first <- combined_data %>%
  filter(format(date, "%m-%d") != "02-29") %>%  
  mutate(
    doy_actual = yday(date),  
    leap_year = lubridate::leap_year(date), 
    Time = ifelse(leap_year & doy_actual > 59, doy_actual - 1, doy_actual) 
  ) %>%
  rename(GPP_mod = GPP) %>%
  group_by(site, Time) %>%
  summarize(
    mean_gpp = mean(GPP_mod, na.rm = T),
    sd_gpp = sd(GPP_mod, na.rm = T)
  ) %>%
  ungroup()




flux_eos_mean_first_then_min2 <- flux_eos_mean_first_then_min %>%
  mutate(
    temp_eos = as.Date(paste(2001, mean_eos, sep = "-"), format = "%Y-%j"),
    temp_sos = as.Date(paste(2001, mean_sos, sep = "-"), format = "%Y-%j"),
    adjusted_sos = as.Date(ifelse(mean_sos < mean_eos, as.Date(paste(2002, mean_sos, sep = "-"), format = "%Y-%j"), temp_sos))
  ) %>%
  filter(!is.na(mean_eos) & !is.na(mean_sos))


doy_sequences_flux <- flux_eos_mean_first_then_min2 %>%
  filter(!is.na(mean_sos) & !is.na(mean_eos)) %>%
  rowwise() %>%
  mutate(
    mean_sos = as.integer(mean_sos),
    mean_eos = as.integer(mean_eos),
    doy_range = list(
      if (mean_sos > mean_eos) {
        c(seq(from = mean_sos, to = 365, by = 1), seq(from = 1, to = mean_eos, by = 1))
      } else {
        seq(from = mean_sos, to = mean_eos, by = 1)
      }
    )
  ) %>%
  ungroup()

doy_sequences_expanded_flux <- doy_sequences_flux %>%
  select(site, doy_range) %>%
  unnest(doy_range) %>%
  mutate(site = str_replace(site, "-", "_")) 


filtered_gpp_flux_eos <- flux_eos_gpp_mean_first %>%
  anti_join(doy_sequences_expanded_flux, by = c("site", "Time" = "doy_range"))


temp_flux_eos <- flux_eos_mean_first_then_min2 %>%
  select(site, mean_eos, mean_sos) %>%
  mutate(mean_eos = as.integer(mean_eos),
         mean_sos = as.integer(mean_sos),
         site = str_replace(site, "-", "_"),
         eos_year = ifelse(mean_sos < mean_eos, 2001, 2001),
         sos_year = ifelse(mean_sos < mean_eos, 2002, 2001),
         mean_eos = paste(mean_eos, eos_year, sep = "-"),
         mean_sos = paste(mean_sos, sos_year, sep = "-")
  ) %>%
  select(-eos_year, -sos_year)  



temp_flux_eos <- temp_flux_eos %>%
  separate(mean_eos, into = c("eos_doy", "eos_year"), sep = "-", convert = TRUE) %>%
  separate(mean_sos, into = c("sos_doy", "sos_year"), sep = "-", convert = TRUE)

expanded_dates_flux <- temp_flux_eos %>%
  mutate(data = pmap(list(site, eos_doy, eos_year, sos_doy, sos_year), 
                     ~ data.frame(
                       site = ..1,
                       doy = if(..5 > ..3) c(seq(..2, 365), seq(1, ..4)) else seq(..2, ..4),
                       year = if(..5 > ..3) c(rep(..3, 365 - ..2 + 1), rep(..5, ..4)) else rep(..3, ..4 - ..2 + 1)
                     ))) %>%
  select(data) %>%
  unnest(data)

filtered_gpp_flux_eos_years <- filtered_gpp_flux_eos %>%
  left_join(expanded_dates_flux, by = c("site", "Time" = "doy")) %>%
  mutate(date = as.Date(paste(year, Time, sep = "-"), format = "%Y-%j")) %>%
  arrange(site, date)


################################################################
# prep and extract senescence period from dynamic model
################################################################


new_eos_mean_first_then_min <- grass_new_eos_long %>% 
  anti_join(invalid_years, by = c("site", "origin")) %>%
  filter(!is.na(value)) %>%
  mutate(value = ifelse(value < 0, 365 + value, value),  
         value = ifelse(value > 365, value - 365, value)) %>% 
  pivot_wider(names_from = metric, values_from = value) %>%
  group_by(site) %>%
  summarize(mean_eos = mean(DER.eos, na.rm =T),
            sd_eos = sd(DER.eos, na.rm =T),
            mean_sos = mean(DER.sos, na.rm =T),
            sd_sos = sd(DER.sos, na.rm =T)) %>%
  ungroup()

new_eos_gpp_mean_first <- matched_data %>%
  select(c("GPP.x", "Site", "Time")) %>%
  rename(GPP_mod = GPP.x, 
         site = Site) %>%
  group_by(site, Time) %>%
  summarize(mean_gpp = mean(GPP_mod, na.rm =T),
            sd_gpp = sd(GPP_mod, na.rm =T)) %>%
  ungroup()

new_eos_mean_first_then_min2 <- new_eos_mean_first_then_min %>%
  mutate(
    temp_eos = as.Date(paste(2001, mean_eos, sep = "-"), format = "%Y-%j"),
    temp_sos = as.Date(paste(2001, mean_sos, sep = "-"), format = "%Y-%j"),
    adjusted_sos = as.Date(ifelse(mean_sos < mean_eos, as.Date(paste(2002, mean_sos, sep = "-"), format = "%Y-%j"), temp_sos))
  ) %>%
  filter(!is.na(mean_eos) & !is.na(mean_sos))


doy_sequences_new <- new_eos_mean_first_then_min2 %>%
  filter(!is.na(mean_sos) & !is.na(mean_eos)) %>%
  rowwise() %>%
  mutate(
    mean_sos = as.integer(mean_sos),
    mean_eos = as.integer(mean_eos),
    doy_range = list(
      if (mean_sos > mean_eos) {
        c(seq(from = mean_sos, to = 365, by = 1), seq(from = 1, to = mean_eos, by = 1))
      } else {
        seq(from = mean_sos, to = mean_eos, by = 1)
      }
    )
  ) %>%
  ungroup()

doy_sequences_expanded_new <- doy_sequences_new %>%
  select(site, doy_range) %>%
  unnest(doy_range) %>%
  mutate(site = str_replace(site, "-", "_")) 

filtered_gpp_new_eos <- new_eos_gpp_mean_first %>%
  anti_join(doy_sequences_expanded_new, by = c("site", "Time" = "doy_range"))

temp_new_eos <- new_eos_mean_first_then_min2 %>%
  select(site, mean_eos, mean_sos) %>%
  mutate(mean_eos = as.integer(mean_eos),
         mean_sos = as.integer(mean_sos),
         site = str_replace(site, "-", "_"),
         eos_year = ifelse(mean_sos < mean_eos, 2001, 2001),
         sos_year = ifelse(mean_sos < mean_eos, 2002, 2001),
         mean_eos = paste(mean_eos, eos_year, sep = "-"),
         mean_sos = paste(mean_sos, sos_year, sep = "-")
  ) %>%
  select(-eos_year, -sos_year) 

temp_new_eos <- temp_new_eos %>%
  separate(mean_eos, into = c("eos_doy", "eos_year"), sep = "-", convert = TRUE) %>%
  separate(mean_sos, into = c("sos_doy", "sos_year"), sep = "-", convert = TRUE)

expanded_dates_new <- temp_new_eos %>%
  mutate(data = pmap(list(site, eos_doy, eos_year, sos_doy, sos_year), 
                     ~ data.frame(
                       site = ..1,
                       doy = if(..5 > ..3) c(seq(..2, 365), seq(1, ..4)) else seq(..2, ..4),
                       year = if(..5 > ..3) c(rep(..3, 365 - ..2 + 1), rep(..5, ..4)) else rep(..3, ..4 - ..2 + 1)
                     ))) %>%
  select(data) %>%
  unnest(data)

filtered_gpp_new_eos_years <- filtered_gpp_new_eos %>%
  left_join(expanded_dates_new, by = c("site", "Time" = "doy")) %>%
  mutate(date = as.Date(paste(year, Time, sep = "-"), format = "%Y-%j")) %>%
  arrange(site, date)


min_gpp_dates_new <- filtered_gpp_new_eos_years %>%
  arrange(site, date) %>% 
  group_by(site) %>% 
  mutate(ma_gpp = rollmean(mean_gpp, 14, fill = NA, align = "center")) %>%  
  ungroup() %>%
  group_by(site) %>%
  filter(ma_gpp == min(ma_gpp, na.rm = TRUE)) %>%
  slice(1) %>%  
  select(site, date, ma_gpp)

filtered_gpp_new_eos_years_after_min <- filtered_gpp_new_eos_years %>%
  inner_join(min_gpp_dates_new %>% select(site, date), by = "site")  %>%
  filter(date.x <= date.y) %>%
  select(-date.y) %>%
  rename(date=date.x)

min_gpp_dates_flux <- filtered_gpp_new_eos_years %>%
  arrange(site, date) %>%  
  group_by(site) %>%  
  mutate(ma_gpp = rollmean(mean_gpp, 14, fill = NA, align = "center")) %>%  
  ungroup() %>%
  group_by(site) %>%
  filter(ma_gpp == min(ma_gpp, na.rm = TRUE)) %>%
  select(site, date, ma_gpp)


filtered_gpp_flux_eos_years_after_min <- filtered_gpp_flux_eos_years %>%
  inner_join(min_gpp_dates_flux %>% select(site, date), by = "site")  %>%
  filter(date.x <= date.y) %>%
  select(-date.y) %>%
  rename(date=date.x)

################################################################
# prep and extract senescence period from default model
################################################################


dev_eos_mean_first_then_min <- grass_dev_eos_long %>% 
  filter(!is.na(value)) %>%
  anti_join(invalid_years_dev, by = c("site", "origin")) %>%
  mutate(value = ifelse(value < 0, 365 + value, value),  
         value = ifelse(value > 365, value - 365, value)) %>% 
  pivot_wider(names_from = metric, values_from = value) %>%
  group_by(site) %>%
  summarize(mean_eos = mean(DER.eos, na.rm =T),
            sd_eos = sd(DER.eos, na.rm =T),
            mean_sos = mean(DER.sos, na.rm =T),
            sd_sos = sd(DER.sos, na.rm =T)) %>%
  ungroup()

dev_eos_gpp_mean_first <- matched_data_dev %>%
  select(c("GPP.x", "Site", "Time")) %>%
  rename(GPP_mod = GPP.x, 
         site = Site) %>%
  group_by(site, Time) %>%
  summarize(mean_gpp = mean(GPP_mod, na.rm =T),
            sd_gpp = sd(GPP_mod, na.rm =T)) %>%
  ungroup()

dev_eos_mean_first_then_min2 <- dev_eos_mean_first_then_min %>%
  mutate(
    temp_eos = as.Date(paste(2001, mean_eos, sep = "-"), format = "%Y-%j"),
    temp_sos = as.Date(paste(2001, mean_sos, sep = "-"), format = "%Y-%j"),
    adjusted_sos = as.Date(ifelse(mean_sos < mean_eos, as.Date(paste(2002, mean_sos, sep = "-"), format = "%Y-%j"), temp_sos))
  ) %>%
  filter(!is.na(mean_eos) & !is.na(mean_sos))

doy_sequences_dev <- dev_eos_mean_first_then_min2 %>%
  filter(!is.na(mean_sos) & !is.na(mean_eos)) %>%
  rowwise() %>%
  mutate(
    mean_sos = as.integer(mean_sos),
    mean_eos = as.integer(mean_eos),
    doy_range = list(
      if (mean_sos > mean_eos) {
        c(seq(from = mean_sos, to = 365, by = 1), seq(from = 1, to = mean_eos, by = 1))
      } else {
        seq(from = mean_sos, to = mean_eos, by = 1)
      }
    )
  ) %>%
  ungroup()

doy_sequences_expanded_dev <- doy_sequences_dev %>%
  select(site, doy_range) %>%
  unnest(doy_range) %>%
  mutate(site = str_replace(site, "-", "_")) 


filtered_gpp_dev_eos <- dev_eos_gpp_mean_first %>%
  anti_join(doy_sequences_expanded_dev, by = c("site", "Time" = "doy_range"))


temp_dev_eos <- dev_eos_mean_first_then_min2 %>%
  select(site, mean_eos, mean_sos) %>%
  mutate(mean_eos = as.integer(mean_eos),
         mean_sos = as.integer(mean_sos),
         site = str_replace(site, "-", "_"),
         eos_year = ifelse(mean_sos < mean_eos, 2001, 2001),
         sos_year = ifelse(mean_sos < mean_eos, 2002, 2001),
         mean_eos = paste(mean_eos, eos_year, sep = "-"),
         mean_sos = paste(mean_sos, sos_year, sep = "-")
  ) %>%
  select(-eos_year, -sos_year) 

temp_dev_eos <- temp_dev_eos %>%
  separate(mean_eos, into = c("eos_doy", "eos_year"), sep = "-", convert = TRUE) %>%
  separate(mean_sos, into = c("sos_doy", "sos_year"), sep = "-", convert = TRUE)

expanded_dates_dev <- temp_dev_eos %>%
  mutate(data = pmap(list(site, eos_doy, eos_year, sos_doy, sos_year), 
                     ~ data.frame(
                       site = ..1,
                       doy = if(..5 > ..3) c(seq(..2, 365), seq(1, ..4)) else seq(..2, ..4),
                       year = if(..5 > ..3) c(rep(..3, 365 - ..2 + 1), rep(..5, ..4)) else rep(..3, ..4 - ..2 + 1)
                     ))) %>%
  select(data) %>%
  unnest(data)

filtered_gpp_dev_eos_years <- filtered_gpp_dev_eos %>%
  left_join(expanded_dates_dev, by = c("site", "Time" = "doy")) %>%
  mutate(date = as.Date(paste(year, Time, sep = "-"), format = "%Y-%j")) %>%
  arrange(site, date)



min_gpp_dates_dev <- filtered_gpp_dev_eos_years %>%
  arrange(site, date) %>%  
  group_by(site) %>%  
  mutate(ma_gpp = rollmean(mean_gpp, 14, fill = NA, align = "center")) %>%
  ungroup() %>%
  group_by(site) %>%
  filter(ma_gpp == min(ma_gpp, na.rm = TRUE)) %>%
  slice(1) %>%  
  select(site, date, ma_gpp)


filtered_gpp_dev_eos_years_after_min <- filtered_gpp_dev_eos_years %>%
  inner_join(min_gpp_dates_dev %>% select(site, date), by = "site")  %>%
  filter(date.x <= date.y) %>%
  select(-date.y) %>%
  rename(date=date.x)

################################################################
# merge all three senescence periods together
################################################################

full_eos_periods_for_all_three <- filtered_gpp_flux_eos_years_after_min %>%
  rename(mean_gpp_flux = mean_gpp, sd_gpp_flux = sd_gpp) %>%  
  inner_join(
    filtered_gpp_new_eos_years_after_min %>%
      rename(mean_gpp_new = mean_gpp, sd_gpp_new = sd_gpp),  
    by = c("site", "date")  
  ) %>%
  inner_join(
    filtered_gpp_dev_eos_years_after_min %>%
      rename(mean_gpp_dev = mean_gpp, sd_gpp_dev = sd_gpp),  
    by = c("site", "date")
  ) %>%
  # REMOVE site AU_DaS, shaky and eos in one week
  filter(site != "AU_DaS")


full_eos_periods_for_dev_and_flux <- filtered_gpp_flux_eos_years_after_min %>%
  rename(mean_gpp_flux = mean_gpp, sd_gpp_flux = sd_gpp) %>%  
  inner_join(
    filtered_gpp_dev_eos_years_after_min %>%
      rename(mean_gpp_dev = mean_gpp, sd_gpp_dev = sd_gpp),  
    by = c("site", "date") 
  )%>%
  # REMOVE site AU_DaS, shaky and eos in one week
  filter(site != "AU_DaS")

full_eos_periods_for_new_and_flux <- filtered_gpp_flux_eos_years_after_min %>%
  rename(mean_gpp_flux = mean_gpp, sd_gpp_flux = sd_gpp) %>% 
  inner_join(
    filtered_gpp_new_eos_years_after_min %>%
      rename(mean_gpp_new = mean_gpp, sd_gpp_new = sd_gpp), 
    by = c("site", "date") 
  )%>%
  filter(site != "AU_DaS")


###########################################################################################################
######################################## CREATE STATS FOR TABLE 6 #########################################
###########################################################################################################

full_eos_periods_for_new_and_flux_gramm <- full_eos_periods_for_new_and_flux %>%
  mutate(across(matches("^mean_gpp_|^sd_gpp_"), ~ . * 12.0107))

full_eos_periods_for_dev_and_flux_gramm <- full_eos_periods_for_dev_and_flux %>%
  mutate(across(matches("^mean_gpp_|^sd_gpp_"), ~ . * 12.0107))

full_eos_periods_for_all_three_gramm <- full_eos_periods_for_all_three %>%
  mutate(across(matches("^mean_gpp_|^sd_gpp_"), ~ . * 12.0107))


dfeos_end <- data.frame(
  metric = c("mae_new", "rmse_new", "nrmse_new", "mae_dev", "rmse_dev", "nrmse_dev"),
  value = c(
    hydroGOF::mae  (obs = full_eos_periods_for_new_and_flux_gramm$mean_gpp_flux, sim = full_eos_periods_for_new_and_flux_gramm$mean_gpp_new),
    hydroGOF::rmse (obs = full_eos_periods_for_new_and_flux_gramm$mean_gpp_flux, sim = full_eos_periods_for_new_and_flux_gramm$mean_gpp_new),
    hydroGOF::nrmse(obs = full_eos_periods_for_new_and_flux_gramm$mean_gpp_flux, sim = full_eos_periods_for_new_and_flux_gramm$mean_gpp_new, norm = 'maxmin'),
    hydroGOF::mae  (obs = full_eos_periods_for_dev_and_flux_gramm$mean_gpp_flux, sim = full_eos_periods_for_dev_and_flux_gramm$mean_gpp_dev),
    hydroGOF::rmse (obs = full_eos_periods_for_dev_and_flux_gramm$mean_gpp_flux, sim = full_eos_periods_for_dev_and_flux_gramm$mean_gpp_dev),
    hydroGOF::nrmse(obs = full_eos_periods_for_dev_and_flux_gramm$mean_gpp_flux, sim = full_eos_periods_for_dev_and_flux_gramm$mean_gpp_dev, norm = 'maxmin')
  )
) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(better.mae = mae_new < mae_dev,
         better.rmse = rmse_new < rmse_dev,
         better.nrmse = nrmse_new < nrmse_dev,
         site = "All Sites") %>%
  select(site, everything())

###########################################################################################################
################################## INDIVIDUAL STATS FOR END PERIOD ########################################
###########################################################################################################

rmse_by_site_new_and_flux <- full_eos_periods_for_new_and_flux_gramm  %>%
  group_by(site) %>%  
  summarise(MAE   = hydroGOF::mae(obs = mean_gpp_flux, sim = mean_gpp_new),
            RMSE  = hydroGOF::rmse(obs = mean_gpp_flux, sim = mean_gpp_new),
            NRMSE = hydroGOF::nrmse(obs = mean_gpp_flux, sim = mean_gpp_new, norm = 'maxmin'))  


rmse_by_site_dev_and_flux <- full_eos_periods_for_dev_and_flux_gramm %>%
  group_by(site) %>% 
  summarise(MAE   = hydroGOF::mae(obs = mean_gpp_flux, sim = mean_gpp_dev),
            RMSE  = hydroGOF::rmse(obs = mean_gpp_flux, sim = mean_gpp_dev),
            NRMSE = hydroGOF::nrmse(obs = mean_gpp_flux, sim = mean_gpp_dev, norm = 'maxmin'))  


rmse_by_site_dev_and_flux_and_new <- rmse_by_site_dev_and_flux  %>%
  inner_join(rmse_by_site_new_and_flux, by = 'site') %>%
  rename(mae_dev = MAE.x,
         mae_new = MAE.y,
         rmse_dev = RMSE.x,
         rmse_new = RMSE.y,
         nrmse_dev = NRMSE.x,
         nrmse_new = NRMSE.y) %>%
  mutate(better.mae = mae_new < mae_dev,
         better.rmse = rmse_new < rmse_dev,
         better.nrmse = nrmse_new < nrmse_dev)

rmse_by_site_dev_and_flux_and_new_all_included <- rmse_by_site_dev_and_flux_and_new %>%
  bind_rows(dfeos_end, .)


dates_for_paper <- full_eos_periods_for_all_three %>%
  select(c(site, date)) %>%
  group_by(site) %>%
  summarize(min_date = min(range(date)),
            max_date = max(range(date))) %>%
  mutate(min_date = yday(min_date),
         max_date = yday(max_date)) %>%
  rename(eos = min_date,
         min_gpp = max_date)


#######################################################################################################
########################################## CREATE FIGURE 6 ############################################
#######################################################################################################


full_eos_periods_for_all_three_temp_plot <- full_eos_periods_for_all_three %>%
  mutate(site = gsub('_' , '-', site))

end_period_all <- ggplot(data = full_eos_periods_for_all_three_temp_plot %>%
                           filter(!(site == "FR-Lq2" & mean_gpp_flux > 0.4162955),
                                  !(site == "US-FPe" & mean_gpp_flux > 0.3))) +
  # Plot the mean GPP flux line and SD ribbon for flux
  geom_line(aes(x = date, y = pmax(0, mean_gpp_flux * 12.0107), color = "EC"), lwd = .5) +
  geom_ribbon(aes(
    x = date,
    ymin = pmax(0, mean_gpp_flux * 12.0107 - sd_gpp_flux * 12.0107),
    ymax = pmax(0, mean_gpp_flux * 12.0107 + sd_gpp_flux * 12.0107)
  ), fill = cbbPalette[4], alpha = 0.3) +
  
  # Plot the mean GPP dev line and SD ribbon for dev
  geom_line(aes(x = date, y = mean_gpp_dev* 12.0107, color = "default"), lwd =.5) +
  geom_ribbon(aes(x = date, ymin = mean_gpp_dev* 12.0107 - sd_gpp_dev* 12.0107, ymax = mean_gpp_dev * 12.0107+ sd_gpp_dev* 12.0107), 
              fill = cbbPalette[6], alpha = 0.3) +  # Add SD ribbon for dev
  
  # Plot the mean GPP dynamic line and SD ribbon for dynamic
  geom_line(aes(x = date, y = mean_gpp_new* 12.0107, color = "dynamic"), lwd =.5) +
  geom_ribbon(aes(x = date, ymin = mean_gpp_new * 12.0107- sd_gpp_new* 12.0107, ymax = mean_gpp_new * 12.0107+ sd_gpp_new* 12.0107), 
              fill = "#ff641e", alpha = 0.3) +  # Add SD ribbon for new
  
  # Facet by site with free axes
  facet_wrap(~site, scales = "free", ncol = 4) +
  # Define color and fill scales
  scale_color_manual(
    name = "Data",
    breaks = c("default", "dynamic", "EC"),
    values = c("default" = cbbPalette[6], "dynamic" = "#ff641e", "EC" = cbbPalette[4]),
    guide = guide_legend(override.aes = list(linewidth = 2))
  ) +
  scale_x_date(
    labels = function(x) str_remove(format(x, "%j"), "^0+"),  # remove leading zeros
    breaks = pretty_breaks(n = 4),
    expand = c(0, 0)
  ) +
  labs(
    title = "",
    x = "",
    y = "GPP g C m<sup>-2</sup> day<sup>-1</sup>") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_text(size = 26, face = 'bold'),
    axis.title.y = element_markdown(size = 26, face = 'plain'),
    axis.title.x = element_text(size = 26, face = 'plain'),
    axis.text.x = element_text(size = 22, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 26, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 26, face = 'plain'),
    legend.title = element_text(size = 26, face = 'bold'),
    strip.text.x = element_text(size = 26, face = 'plain')
  )  ; print(end_period_all)






