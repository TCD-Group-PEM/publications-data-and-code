
 #######################################################################################
 ###### !!! THIS FILE REQUIRES THE DATA LOADED IN FILE 'Figure 3 and Table 4' !!! ######
 #######################################################################################

 
 
######################################################################## 
#################### Create the results for table 2 #################### 
########################################################################
 
################################################ 
#################### US-Tol #################### 
################################################
 
temp_ustol_gpp_sum <- US_ICs_EC %>%
  mutate(Date   = as.Date(Datetime), 
         GPP_NT = GPP_NT * 1800/1e6) %>% # convert from umol s-1 every 30 mins to Mol day
  group_by(Date) %>%
  summarize(GPP_NT = sum(GPP_NT)) %>%
  na.omit()%>%  
  filter(!(format(Date, "%m") == "02" & format(Date, "%d") == "29")) %>%
  ungroup() %>%
  mutate(year   = year(Date)) %>%
  group_by(year) %>%
  summarize(gpp_sum = sum(GPP_NT*12.0107)) %>%
  ungroup() 

temp_ustol_gpp_sum_sim <- USTolsimnormtau_IC_match %>%
  mutate(year   = year(Date)) %>%
  group_by(year) %>%
  summarize(gpp_sum = sum(GPP*12.0107)) %>%
  ungroup() 

temp_ustol_gpp_sum_dev <- USToldev_IC_match %>%
  mutate(year   = year(Date)) %>%
  group_by(year) %>%
  summarize(gpp_sum = sum(GPP*12.0107)) %>%
  ungroup() 

temp_ustol_gpp_sum_dev$Dataset <- "default"
temp_ustol_gpp_sum_sim$Dataset <- "new"
temp_ustol_gpp_sum$Dataset <- "EC"

# Bind them into one dataframe
temp_us_gpp_all_three <- bind_rows(temp_ustol_gpp_sum_dev, 
                                   temp_ustol_gpp_sum_sim, 
                                   temp_ustol_gpp_sum)


################################################ 
#################### ES-LMa #################### 
################################################

temp_esmag_gpp_sum <- maj_ec_flux_total %>%
  mutate(Date   = as.Date(date)) %>% # convert from umol s-1 every 30 mins to Mol day
  group_by(Date) %>%
  summarize(GPP = sum(GPP)) %>%
  na.omit()%>%  
  filter(!(format(Date, "%m") == "02" & format(Date, "%d") == "29")) %>%
  ungroup() %>%
  mutate(year   = year(Date)) %>%
  group_by(year) %>%
  summarize(gpp_sum = sum(GPP*12.0107)) %>%
  ungroup() 

temp_esmag_gpp_sum_sim <- ESMaGsimnormtau %>%
  filter(as.Date(Date) %in% as.Date(maj_ec_flux_total$date)) %>%
  mutate(year   = year(Date)) %>%
  group_by(year) %>%
  summarize(gpp_sum = sum(GPP*12.0107)) %>%
  ungroup() 

temp_esmag_gpp_sum_dev <- ESMaGdev %>%
  filter(as.Date(Date) %in% as.Date(maj_ec_flux_total$date)) %>%
  mutate(year   = year(Date)) %>%
  group_by(year) %>%
  summarize(gpp_sum = sum(GPP*12.0107)) %>%
  ungroup() 

temp_esmag_gpp_sum_dev$Dataset <- "default"
temp_esmag_gpp_sum_sim$Dataset <- "new"
temp_esmag_gpp_sum$Dataset <- "EC"

# Bind them into one dataframe
temp_es_gpp_all_three <- bind_rows(temp_esmag_gpp_sum_dev, 
                                   temp_esmag_gpp_sum_sim, 
                                   temp_esmag_gpp_sum)

temp_es_gpp_all_three <- temp_es_gpp_all_three %>% 
  filter(year > 2003)


################################################ 
#################### IE-Dri #################### 
################################################

temp_iedri_gpp_sum <- plumfluxIEDri %>%
  mutate(Date   = as.Date(date)) %>% # convert from umol s-1 every 30 mins to Mol day
  group_by(Date) %>%
  summarize(GPP = sum(GPP)) %>%
  na.omit()%>%  
  filter(!(format(Date, "%m") == "02" & format(Date, "%d") == "29")) %>%
  ungroup() %>%
  mutate(year   = year(Date)) %>%
  group_by(year) %>%
  summarize(gpp_sum = sum(GPP*12.0107)) %>%
  ungroup() 

temp_iedri_gpp_sum_sim <- IEDrisimnormtau %>%
  filter(as.Date(Date) %in% as.Date(plumfluxIEDri$date)) %>%
  mutate(year   = year(Date)) %>%
  group_by(year) %>%
  summarize(gpp_sum = sum(GPP*12.0107)) %>%
  ungroup() 

temp_iedri_gpp_sum_dev <- IEDridev %>%
  filter(as.Date(Date) %in% as.Date(plumfluxIEDri$date)) %>%
  mutate(year   = year(Date)) %>%
  group_by(year) %>%
  summarize(gpp_sum = sum(GPP*12.0107)) %>%
  ungroup() 

temp_iedri_gpp_sum_dev$Dataset <- "default"
temp_iedri_gpp_sum_sim$Dataset <- "new"
temp_iedri_gpp_sum$Dataset <- "EC"

# Bind them into one dataframe
temp_ie_gpp_all_three <- bind_rows(temp_iedri_gpp_sum_dev, 
                                   temp_iedri_gpp_sum_sim, 
                                   temp_iedri_gpp_sum)


temp_ie_gpp_all_three <- temp_ie_gpp_all_three %>% 
  filter(year > 2002)


temp_ie_gpp_all_three <- temp_ie_gpp_all_three %>%
  mutate(site = "IE-Dri")

temp_es_gpp_all_three <- temp_es_gpp_all_three %>%
  mutate(site = "ES-LMa") 

temp_us_gpp_all_three <- temp_us_gpp_all_three %>%
  mutate(site = "US-Tol")

test_annual_gpp_diff_ie <- aov(data = temp_ie_gpp_all_three, formula = gpp_sum ~ Dataset)
summary.lm(test_annual_gpp_diff_ie)

TukeyHSD(test_annual_gpp_diff_ie)

test_annual_gpp_diff_es <- aov(data = temp_es_gpp_all_three, formula = gpp_sum ~ Dataset)
summary.lm(test_annual_gpp_diff_es)

TukeyHSD(test_annual_gpp_diff_es)

test_annual_gpp_diff_us <- aov(data = temp_us_gpp_all_three, formula = gpp_sum ~ Dataset)
summary.lm(test_annual_gpp_diff_us)

TukeyHSD(test_annual_gpp_diff_us)



########################################################################
######################### Create the Figure 2 ########################## 
########################################################################

three_sites_gpp_for_barplot <- bind_rows(temp_ie_gpp_all_three, temp_es_gpp_all_three, temp_us_gpp_all_three %>% mutate(year = as.numeric((year)))) %>%
  mutate(year = as.factor(year))

mean_gpp_three_sites_values <- three_sites_gpp_for_barplot %>%
  group_by(site, Dataset) %>%
  summarize(gpp_mean = mean(gpp_sum),
            gpp_sd   = sd(gpp_sum),
            .groups = 'drop')


site_results <- data.frame(
  site=c('ES-LMa','ES-LMa','ES-LMa', 'IE-Dri','IE-Dri', 'IE-Dri', 'US-Tol', 'US-Tol', 'US-Tol'),
  data=c('default', 'new', 'EC', 'default', 'new', 'EC', 'default', 'new', 'EC'),
  site_letters=c('a', 'a', 'a', 'a', 'b', 'b', 'a', 'a', 'a')
)

boxplot_stats <- three_sites_gpp_for_barplot %>%
  group_by(site, Dataset) %>%
  summarise(
    y_max = max(gpp_sum) + 10,
    n = n(),
    .groups = 'drop'
  ) %>%
  left_join(site_results, by = c("site", "Dataset"='data'))


# Add the significance letters and positions
three_sites_gpp_for_barplot2 <- three_sites_gpp_for_barplot %>%
  left_join(boxplot_stats, by = c("site", "Dataset")) %>%
  mutate(Dataset = factor(Dataset, levels = c("EC", "default", "new")))


# Create the ggplot with significance letters
three_sites_annual_gpp_figure <- ggplot(three_sites_gpp_for_barplot2, aes(x = Dataset, y = gpp_sum, fill = Dataset)) +
  geom_boxplot(position = 'dodge', width = 0.7) +
  facet_wrap(~site, scales = 'free') +  # Create separate plots for each site
  labs(x = "", y = "Annual GPP (g C m<sup>-2</sup> year<sup>-1</sup>)") + 
  scale_fill_manual(
    name = "Data", 
    breaks = c("default", "new", "EC"), 
    labels = c("default" = "default", "new" = "dynamic", "EC" = "EC"),
    values = c("default" = cbbPalette[6], "new" = "#ff641e", "EC" = cbbPalette[4])
  ) + 
  scale_x_discrete(labels = c("default" = "default", "new" = "dynamic", "EC" = "EC")) +
  geom_text(aes(label = site_letters, y = y_max), 
            position = position_dodge(width = 0.7), 
            vjust = -0.55,
            size = 12) +  
  geom_text(data = boxplot_stats, aes(x = Inf, y = Inf, label = paste("n =", n)), 
            hjust = 1.2, vjust = 1.2, size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.title.y = element_markdown(size = 36),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold'),
    strip.text.x = element_text(size = 36, face = 'bold')
  ); print(three_sites_annual_gpp_figure)

