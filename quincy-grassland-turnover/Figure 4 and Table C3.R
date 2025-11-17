
library(tidyverse)
library(ggtext)
library(zoo)
library(cowplot)
library(grid)
library(ggpubr)
library(purrr)

# create color blind friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# create helper day of the year vector
yearvec12spin<-rep(c(1:150),each=365)


#################################################################################################
######################################### LOAD THE DATA #########################################
#################################################################################################

ecoC_ESLMa_new <- read.delim("/Users/josuaseitz/remote/PEM/users/seitzj/Chapter_1/simulations/eco_c_turnover/ES-LMa/newT_ecoC.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "ES-LMa_new" )

ecoC_ESLMa_dev <- read.delim("/Users/josuaseitz/remote/PEM/users/seitzj/Chapter_1/simulations/eco_c_turnover/ES-LMa/dev_ecoC.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "ES-LMa_dev" )

ecoC_IEDri_new <- read.delim("/Users/josuaseitz/remote/PEM/users/seitzj/Chapter_1/simulations/eco_c_turnover/IE-Dri/newT_ecoC.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "IE-Dri_new" )

ecoC_IEDri_dev <- read.delim("/Users/josuaseitz/remote/PEM/users/seitzj/Chapter_1/simulations/eco_c_turnover/IE-Dri/dev_ecoC.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "IE-Dri_dev" )

ecoC_USTOL_new <- read.delim("/Users/josuaseitz/remote/PEM/users/seitzj/Chapter_1/simulations/eco_c_turnover/US-TOL/newT_ecoC.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "US-Tol_new" )

ecoC_USTOL_dev <- read.delim("/Users/josuaseitz/remote/PEM/users/seitzj/Chapter_1/simulations/eco_c_turnover/US-TOL/dev_ecoC.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "US-Tol_dev" )

ecoC_CGE_new <- read.delim("/Users/josuaseitz/remote/PEM/users/seitzj/Chapter_1/simulations/eco_c_turnover/CGE/newT_ecoC.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "CGE_new" )

ecoC_CGE_dev <- read.delim("/Users/josuaseitz/remote/PEM/users/seitzj/Chapter_1/simulations/eco_c_turnover/CGE/dev_ecoC.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "CGE_dev" )


EcoC_combined <- bind_rows(
  ecoC_CGE_dev %>%
    mutate(site = 'CGE',
           Source = 'default model'),
  ecoC_CGE_new %>%
    mutate(site = 'CGE',
           Source = 'new model'),
  ecoC_USTOL_dev %>%
    mutate(site = 'US-Tol',
           Source = 'default model'),
  ecoC_USTOL_new %>%
    mutate(site = 'US-Tol',
           Source = 'new model'),
  ecoC_IEDri_dev %>%
    mutate(site = 'IE-Dri',
           Source = 'default model'),
  ecoC_IEDri_new %>%
    mutate(site = 'IE-Dri',
           Source = 'new model'), 
  ecoC_ESLMa_dev %>%
    mutate(site = 'ES-LMa',
           Source = 'default model'),
  ecoC_ESLMa_new %>%
    mutate(site = 'ES-LMa',
           Source = 'new model')
)


###################################################################################################
######################################### Create Figure 4 #########################################
###################################################################################################


EcoC_combined <- EcoC_combined %>% 
  filter(year %in% 100:122) %>%
  mutate(Source = ifelse(Source == "new model", "dynamic model", Source))

p1 <- ggplot() +
  geom_boxplot(data = EcoC_combined , 
               aes(x=Source, y = EcoC * 12.0107 /1000, fill = Source)) +
  facet_wrap(~site, scales = 'free') +
  labs(
    title = "Ecosystem C at the four main sites 2000-2022",
    x = "",
    y = "Eco C (kg C m<sup>-2</sup>)")+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    text = element_text(size = 14),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown()
  ) + 
  scale_fill_manual(
    name = "Data",
    values = c(
      "default model" = cbbPalette[6], 
      "dynamic model" = "#ff641e")
  ) ; print(p1)


p2 <- ggplot() +
  geom_boxplot(data = EcoC_combined , aes(x=Source, y = (TotalSoilOrgC) * 12.0107 /1000, fill = Source)) +
  facet_wrap(~site, scales = 'free') +
  labs(
    title = "SOC at the four main sites 2000-2022",
    x = "",
    y = "SOC (kg C m<sup>-2</sup>)")+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    text = element_text(size = 14),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown()
  ) + 
  scale_fill_manual(
    name = "Data",
    values = c(
      "default model" = cbbPalette[6], 
      "dynamic model" = "#ff641e")
  ) ; print(p2)



p3 <- ggplot() +
  geom_boxplot(data = EcoC_combined , aes(x=Source, y = (EcoC - TotalSoilOrgC) * 12.0107 /1000, fill = Source)) +
  facet_wrap(~site, scales = 'free') +
  labs(
    title = "Vegetation C at the four main sites 2000-2022",
    x = "",
    y = "Vegetation C (kg C m<sup>-2</sup>)")+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    text = element_text(size = 14),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown()
  ) + 
  scale_fill_manual(
    name = "Data",
    values = c(
      "default model" = cbbPalette[6], 
      "dynamic model" = "#ff641e")
  ) ; print(p3)


p4 <- ggarrange(p1 + theme(legend.position = "none"), p2+ theme(legend.position = "none"), p3+ theme(legend.position = "bottom"), ncol = 1 ) ; p4



EcoC_yearly_mean <- EcoC_combined %>%
  filter(year %in% 100:122) %>%
  group_by(year, site, Source) %>%
  summarize(
    SoilC = mean(TotalSoilOrgC * 12.0107 / 1000, na.rm = TRUE),
    VegC = mean(EcoC * 12.0107 / 1000 - TotalSoilOrgC * 12.0107 / 1000 , na.rm = TRUE),
    EcoC = mean(EcoC * 12.0107 / 1000, na.rm = TRUE),
    VegC_sd = sd(VegC * 12.0107 / 1000 - TotalSoilOrgC * 12.0107 / 1000 , na.rm = TRUE),
    SoilC_sd = sd(TotalSoilOrgC * 12.0107 / 1000, na.rm = TRUE),
    EcoC_sd = sd(EcoC * 12.0107 / 1000 - TotalSoilOrgC * 12.0107 / 1000 , na.rm = TRUE),
    .groups = "drop"
  )%>%
  mutate(Source = recode(Source,
                         "default model" = "default",
                         "dynamic model" = "dynamic",
                         "EC flux" = "EC"))


p5 <- ggplot(EcoC_yearly_mean) +
  geom_boxplot(
    aes(x = Source, y = EcoC, fill = Source),
  ) +
  geom_point(
    aes(x = Source, y = EcoC), color = 'black', alpha = .4) +
  facet_wrap(~site, scales = 'free', ncol = 4) +
  labs(
    title = "mean annual Ecosystem C at the four main sites 2000-2022",
    x = "",
    y = "Ecosystem C (kg C m<sup>-2</sup>)")+
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.x = element_text(size = 26, face = 'plain'),
    axis.title.y = element_markdown(size = 26),
    axis.text.x = element_text(size = 26, face = 'plain'),   
    axis.text.y = element_text(size = 26, face = 'plain'),   
    legend.text = element_text(size = 26, face = 'plain'),
    legend.title = element_text(size = 26, face = 'bold'),
    strip.text.x = element_text(size = 26, face = 'bold')
  )  + 
  scale_fill_manual(
    name = "Data",
    values = c(
      "default" = cbbPalette[6], 
      "dynamic" = "#ff641e")
  ) ; print(p5)

p6 <- ggplot(EcoC_yearly_mean) +
  geom_boxplot(
    aes(x = Source, y = SoilC, fill = Source),
  ) +
  geom_point(
    aes(x = Source, y = SoilC), color = 'black', alpha = .4) +
  facet_wrap(~site, scales = 'free', ncol = 4) +
  labs(
    title = "mean annual SOC at the four main sites 2000-2022",
    x = "",
    y = "SOC (kg C m<sup>-2</sup>)")+
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.x = element_text(size = 26, face = 'plain'),
    axis.title.y = element_markdown(size = 26),
    axis.text.x = element_text(size = 26, face = 'plain'),  
    axis.text.y = element_text(size = 26, face = 'plain'),  
    legend.text = element_text(size = 26, face = 'plain'),
    legend.title = element_text(size = 26, face = 'bold'),
    strip.text.x = element_text(size = 26, face = 'bold')
  )  + 
  scale_fill_manual(
    name = "Data",
    values = c(
      "default" = cbbPalette[6], 
      "dynamic" = "#ff641e")
  ) ; print(p6)


p7 <- ggplot(EcoC_yearly_mean) +
  geom_boxplot(
    aes(x = Source, y = VegC, fill = Source),
  ) +
  geom_point(
    aes(x = Source, y = VegC), color = 'black', alpha = .4) +
  facet_wrap(~site, scales = 'free', ncol = 4) +
  labs(
    title = "mean annual Vegetation C at the four main sites 2000-2022",
    x = "",
    y = "Vegetation C (kg C m<sup>-2</sup>)")+
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.x = element_text(size = 26, face = 'plain'),
    axis.title.y = element_markdown(size = 26),
    axis.text.x = element_text(size = 26, face = 'plain'),   
    axis.text.y = element_text(size = 26, face = 'plain'),   
    legend.text = element_text(size = 26, face = 'plain'),
    legend.title = element_text(size = 26, face = 'bold'),
    strip.text.x = element_text(size = 26, face = 'bold')
  ) + 
  scale_fill_manual(
    name = "Data",
    values = c(
      "default" = cbbPalette[6], 
      "dynamic" = "#ff641e")
  ) ; print(p7)


p8 <- ggarrange( 
  p5 + theme(legend.position = "none"),
  p7 + theme(legend.position = "none"), 
  p6 + theme(legend.position = "none"), nrow =3, 
  labels = c("a)", "b)", 'c)'),  
  label.x = 0,
  label.y = 1,
  font.label = list(size = 36)) ; print(p8)


ggsave(filename = "/path/to/figure.png", 
       plot = p8, 
       scale = 5,
       width = 1000,
       height = 1000,
       units = "px")



###################################################################################################
######################################### Create Table C3 #########################################
###################################################################################################

vars <- c("EcoC", "VegC", "SoilC")


results <- EcoC_yearly_mean %>%
  group_by(site) %>%
  group_split() %>%
  map_df(function(df_site) {
    site_name <- unique(df_site$site)
    
    map_df(vars, function(var) {
      dynamic_vals <- df_site %>% filter(Source == "dynamic") %>% pull(all_of(var))
      default_vals <- df_site %>% filter(Source == "default") %>% pull(all_of(var))
      
      if (length(dynamic_vals) >= 2 && length(default_vals) >= 2) {
        t_res <- t.test(dynamic_vals, default_vals)
        tibble(
          site = site_name,
          variable = var,
          t_statistic = t_res$statistic,
          df = t_res$parameter,
          p_value = t_res$p.value,
          mean_dynamic = mean(dynamic_vals, na.rm = TRUE),
          mean_default = mean(default_vals, na.rm = TRUE)
        )
      } else {
        tibble(
          site = site_name,
          variable = var,
          t_statistic = NA,
          df = NA,
          p_value = NA,
          mean_dynamic = mean(dynamic_vals, na.rm = TRUE),
          mean_default = mean(default_vals, na.rm = TRUE)
        )
      }
    })
  })

print(results)

results_df <- as.data.frame(results) %>%
  mutate(
    p_value = ifelse(p_value < 0.001, "<0.001", p_value)
  ) %>%
  mutate(across(where(is.numeric) & !matches("p_value"), ~ round(.x, 2)))
