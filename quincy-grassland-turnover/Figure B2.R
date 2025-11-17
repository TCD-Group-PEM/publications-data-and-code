
# load libraries
library(ggtext)
library(tidyverse)



########################################################################
################################ SETUP #################################
########################################################################

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



########################################################################
########################### LOAD THE DATA ##############################
########################################################################
ustol_stdout_turn <- read.delim('/path/to/quincy/dynamic/US-Tol/std_out/file.out', header = F, sep = "", skip = 15)

es_stdout_turn <- read.delim('/path/to/quincy/dynamic/ES-LMa/std_out/file.out', header = F, sep = "", skip = 15)

ie_stdout_turn <- read.delim('/path/to/quincy/dynamic/IE-Dri/std_out/file.out', header = F, sep = "", skip = 15)

cg_stdout_turn <- read.delim('/path/to/quincy/dynamic/CGE/std_out/file.out', header = F, sep = "", skip = 15)


##################################
# US-Tol
##################################

usvec<- row.names(ustol_stdout_turn)

ustol_stdout_turn2 <- ustol_stdout_turn %>% select(c(1,2))

colnames(ustol_stdout_turn2)[2] <- 'val'
colnames(ustol_stdout_turn2)[1] <- 'type'

ustol_stdout_turn3 <- ustol_stdout_turn2 %>%
  mutate(val = as.numeric(val),
         type = as.factor(type)) %>%
  filter(type %in% c('temppart:', 'moistpart:', 'lightpart:', 'shedval:', 'kstar_labile:')) %>%
  mutate(type = droplevels(type),
         timestep = rep(1:(n() / 5), each = 5),
         year = if_else(timestep <= 365 * 48, 1, 2),
         day_of_year = rep(rep(1:365, each = 48), times = 2) %>% rep(each = 5),
         time_step_in_day = rep(rep(0:47, times = 365 * 2), each = 5),
         hour = floor(time_step_in_day / 2),
         minute = if_else(time_step_in_day %% 2 == 0, 0, 30),
         date_time = ymd_hm(paste0("200", year, "-01-01 ", sprintf("%02d:%02d", hour, minute))) + days(day_of_year - 1)
  ) %>%
  select(-timestep, -time_step_in_day, -hour, -minute)


ustol_stdout_turn4 <- ustol_stdout_turn3 %>%
  group_by(day_of_year, type) %>%
  summarize(val = mean(val, na.rm = T)) %>%
  ungroup()


uscomp_plot <- ggplot()+
  #geom_line(data=ustol_stdout_turn2 %>% filter(type %in% 'fturn_leaf:'), aes(x=ts, y = val, color = 'turn')) +
  geom_line(data=ustol_stdout_turn4 %>% filter(type %in% 'lightpart:'), aes(x=day_of_year, y = 1-val, color = 'Light'), lwd = 1) +
  geom_line(data=ustol_stdout_turn4 %>% filter(type %in% 'temppart:'), aes(x=day_of_year, y = 1-val, color = 'Temperature'), lwd = 1) +
  geom_line(data=ustol_stdout_turn4 %>% filter(type %in% 'moistpart:'), aes(x=day_of_year, y = 1-val, color = 'Moisture'), lwd = 1) +
  geom_line(data=ustol_stdout_turn4 %>% filter(type %in% 'shedval:'), aes(x=day_of_year, y = val, color = 'Combined'), lwd = 1) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    text = element_text(size = 10),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown()
  ) +  
  scale_color_manual(
    name = "Model Component",
    breaks = c("Light", "Temperature", "Moisture", "Combined"),
    values = c("Light" = cbbPalette[2], "Temperature" = 'red', "Moisture" = cbbPalette[6], "Combined" = cbbPalette[4])
  ) +
  labs(
    title = "d)    US-Tol",
    x = "Day of Year",
    y = "Impact")



##################################
# ES-LMa
##################################

esvec<- row.names(es_stdout_turn)

es_stdout_turn2 <- es_stdout_turn %>% select(c(1,2))

colnames(es_stdout_turn2)[2] <- 'val'
colnames(es_stdout_turn2)[1] <- 'type'

es_stdout_turn3 <- es_stdout_turn2 %>%
  mutate(val = as.numeric(val),
         type = as.factor(type)) %>%
  filter(type %in% c('temppart:', 'moistpart:', 'lightpart:', 'shedval:', 'kstar_labile:')) %>%
  mutate(type = droplevels(type),
         timestep = rep(1:(n() / 5), each = 5),
         year = if_else(timestep <= 365 * 48, 1, 2),
         day_of_year = rep(rep(1:365, each = 48), times = 2) %>% rep(each = 5),
         time_step_in_day = rep(rep(0:47, times = 365 * 2), each = 5),
         hour = floor(time_step_in_day / 2),
         minute = if_else(time_step_in_day %% 2 == 0, 0, 30),
         date_time = ymd_hm(paste0("200", year, "-01-01 ", sprintf("%02d:%02d", hour, minute))) + days(day_of_year - 1)
  ) %>%
  select(-timestep, -time_step_in_day, -hour, -minute)


es_stdout_turn4 <- es_stdout_turn3 %>%
  group_by(day_of_year, type) %>%
  summarize(val = mean(val, na.rm = T)) %>%
  ungroup()


escomp_plot <- ggplot()+
  #geom_line(data=es_stdout_turn2 %>% filter(type %in% 'fturn_leaf:'), aes(x=ts, y = val, color = 'turn')) +
  geom_line(data=es_stdout_turn4 %>% filter(type %in% 'lightpart:'), aes(x=day_of_year, y = 1-val, color = 'Light'), lwd = 1) +
  geom_line(data=es_stdout_turn4 %>% filter(type %in% 'temppart:'), aes(x=day_of_year, y = 1-val, color = 'Temperature'), lwd = 1) +
  geom_line(data=es_stdout_turn4 %>% filter(type %in% 'moistpart:'), aes(x=day_of_year, y = 1-val, color = 'Moisture'), lwd = 1) +
  geom_line(data=es_stdout_turn4 %>% filter(type %in% 'shedval:'), aes(x=day_of_year, y = val, color = 'Combined'), lwd = 1) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    text = element_text(size = 10),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown()
  ) +  
  scale_color_manual(
    name = "Model Component",
    breaks = c("Light", "Temperature", "Moisture", "Combined"),
    values = c("Light" = cbbPalette[2], "Temperature" = 'red', "Moisture" = cbbPalette[6], "Combined" = cbbPalette[4])
  ) +
  labs(
    title = "a)    ES-LMa",
    x = "Day of Year",
    y = "Impact")


##################################
# IE-Dri
##################################

ievec<- row.names(ie_stdout_turn)

ie_stdout_turn2 <- ie_stdout_turn %>% select(c(1,2))

colnames(ie_stdout_turn2)[2] <- 'val'
colnames(ie_stdout_turn2)[1] <- 'type'

ie_stdout_turn3 <- ie_stdout_turn2 %>%
  mutate(val = as.numeric(val),
         type = as.factor(type)) %>%
  filter(type %in% c('temppart:', 'moistpart:', 'lightpart:', 'shedval:', 'kstar_labile:')) %>%
  mutate(type = droplevels(type),
         timestep = rep(1:(n() / 5), each = 5),
         year = if_else(timestep <= 365 * 48, 1, 2),
         day_of_year = rep(rep(1:365, each = 48), times = 2) %>% rep(each = 5),
         time_step_in_day = rep(rep(0:47, times = 365 * 2), each = 5),
         hour = floor(time_step_in_day / 2),
         minute = if_else(time_step_in_day %% 2 == 0, 0, 30),
         date_time = ymd_hm(paste0("200", year, "-01-01 ", sprintf("%02d:%02d", hour, minute))) + days(day_of_year - 1)
  ) %>%
  select(-timestep, -time_step_in_day, -hour, -minute)


ie_stdout_turn4 <- ie_stdout_turn3 %>%
  group_by(day_of_year, type) %>%
  summarize(val = mean(val, na.rm = T)) %>%
  ungroup()


iecomp_plot <- ggplot()+
  geom_line(data=ie_stdout_turn4 %>% filter(type %in% 'lightpart:'), aes(x=day_of_year, y = 1-val, color = 'Light'), lwd = 1) +
  geom_line(data=ie_stdout_turn4 %>% filter(type %in% 'temppart:'), aes(x=day_of_year, y = 1-val, color = 'Temperature'), lwd = 1) +
  geom_line(data=ie_stdout_turn4 %>% filter(type %in% 'moistpart:'), aes(x=day_of_year, y = 1-val, color = 'Moisture'), lwd = 1) +
  geom_line(data=ie_stdout_turn4 %>% filter(type %in% 'shedval:'), aes(x=day_of_year, y = val, color = 'Combined'), lwd = 1) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    text = element_text(size = 10),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown()
  ) +  
  scale_color_manual(
    name = "Model Component",
    breaks = c("Light", "Temperature", "Moisture", "Combined"),
    values = c("Light" = cbbPalette[2], "Temperature" = 'red', "Moisture" = cbbPalette[6], "Combined" = cbbPalette[4])
  ) +
  labs(
    title = "c)    IE-Dri",
    x = "Day of Year",
    y = "Impact")



##################################
# CGE
##################################

cgvec<- row.names(cg_stdout_turn)

cg_stdout_turn2 <- cg_stdout_turn %>% select(c(1,2))

colnames(cg_stdout_turn2)[2] <- 'val'
colnames(cg_stdout_turn2)[1] <- 'type'

cg_stdout_turn3 <- cg_stdout_turn2 %>%
  mutate(val = as.numeric(val),
         type = as.factor(type)) %>%
  filter(type %in% c('temppart:', 'moistpart:', 'lightpart:', 'shedval:', 'kstar_labile:')) %>%
  mutate(type = droplevels(type),
         timestep = rep(1:(n() / 5), each = 5),
         year = if_else(timestep <= 365 * 48, 1, 2),
         day_of_year = rep(rep(1:365, each = 48), times = 2) %>% rep(each = 5),
         time_step_in_day = rep(rep(0:47, times = 365 * 2), each = 5),
         hour = floor(time_step_in_day / 2),
         minute = if_else(time_step_in_day %% 2 == 0, 0, 30),
         date_time = ymd_hm(paste0("200", year, "-01-01 ", sprintf("%02d:%02d", hour, minute))) + days(day_of_year - 1)
  ) %>%
  select(-timestep, -time_step_in_day, -hour, -minute)


cg_stdout_turn4 <- cg_stdout_turn3 %>%
  group_by(day_of_year, type) %>%
  summarize(val = mean(val, na.rm = T)) %>%
  ungroup()

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cgcomp_plot <- ggplot()+
  geom_line(data=cg_stdout_turn4 %>% filter(type %in% 'lightpart:'), aes(x=day_of_year, y = 1-val, color = 'Light'), lwd =1) +
  geom_line(data=cg_stdout_turn4 %>% filter(type %in% 'temppart:'), aes(x=day_of_year, y = 1-val, color = 'Temperature'), lwd =1) +
  geom_line(data=cg_stdout_turn4 %>% filter(type %in% 'moistpart:'), aes(x=day_of_year, y = 1-val, color = 'Moisture'), lwd =1) +
  geom_line(data=cg_stdout_turn4 %>% filter(type %in% 'shedval:'), aes(x=day_of_year, y = val, color = 'Combined'), lwd =1) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    text = element_text(size = 10),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown()
  ) +  
  scale_color_manual(
    name = "Model Component",
    breaks = c("Light", "Temperature", "Moisture", "Combined"),
    values = c("Light" = cbbPalette[2], "Temperature" = 'red', "Moisture" = cbbPalette[6], "Combined" = cbbPalette[4])
  ) +
  labs(
    title = "b)    CGE",
    x = "Day of Year",
    y = "Impact"); print(cgcomp_plot)


allresponses_plot <- ggarrange(escomp_plot, cgcomp_plot, iecomp_plot, uscomp_plot, ncol = 2, nrow = 2 , common.legend = T, legend = 'bottom')
print(allresponses_plot)

