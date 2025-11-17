
# GPP data:
# ES-LMa: 2003-2006 (Source: PLUMBER2), 2015-2018 (Source: Nair et al. 2024)
# US-Tol: 2014-2016 (Source: FLUXNET site US-ICs)  
# IE-Dri: 2002-2005 (Source: PLUMBER2)


 # load libraries
library(tidyverse)
library(ncdf4)
library(ggtext)
library(phenocamapi)
library(hydroGOF) 
library(INDperform)


 # create color blind friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



########################################################################
################################ SETUP #################################
########################################################################

 # PLUMBER2 scripts for obtaining flux and met data: https://github.com/aukkola/PLUMBER2

 # function for processing GPP data in netcdf format to convert into df and add Date column + convert GPP units
process_nc_file <- function(file_path) {
  file_name <- basename(file_path)
  site_name <- strsplit(file_name, "_")[[1]][1]  
  year <- substr(strsplit(file_name, "_")[[1]][2], 1, 4)  
  
  site_name <- gsub("-", "_", site_name)
  
  nc <- nc_open(file_path)
  
  print(nc)
  
  variable_names <- names(nc$var)
  
  GPP <- ncvar_get(nc, "GPP")
  time <- ncvar_get(nc, "time")
  
  nc_close(nc)
  
  origin_date <- as.POSIXct(paste0(year, "-01-01 00:00:00"), tz = "CET")
  
  date <- origin_date + time
  
  data <- data.frame(date = date, GPP = GPP)
  
  df_name <- paste0(site_name, "EC")
  assign(df_name, data, envir = .GlobalEnv)
  
  data$GPP_umol_s <- (data$GPP * 1e6) / 86400
  
  data$date <- as.Date(data$date)
  data <- aggregate(GPP ~ date, data = data, FUN = function(x) sum(x * 1800) / 1e6)
  
  assign(df_name, data, envir = .GlobalEnv)
  
  return(data)
}


#################################################################################################
######################################### LOAD THE DATA #########################################
#################################################################################################


########################################################################
################################ ES-LMa ################################
########################################################################


#####################
######## GPP ########
#####################

 # Load QUINCY simulation results and create Date column:
ESMaGsimnormtau <- read.delim("path/to/quincy/dynamic/ES-MaG/output/vegfluxC_daily.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "sim" )

ESMaGdev <- read.delim("path/to/quincy/default/ES-MaG/output/vegfluxC_daily.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "dev" )


# read ES-LMa data from Tower 2015-2018 (DOI = 10.5281/zenodo.1314194) 
# described in Perez-Priego, O. et al 2017 (DOI = 10.1016/j.agrformet.2017.01.009)

subcan<-read.delim("/path/to/file/ESLMa_SubCanopy.txt",sep="")

#convert Julian Day and year to POSIXct date format
subcan$date<-paste(subcan$doy,subcan$year, sep = "/")  %>%
  as.POSIXct(.,format="%j/%Y")

#select only GPP and Time and transform from half-hourly to daily
flux2<-subcan[,c(23,26)] %>%
  group_by(date)%>%
  summarise(GPP = sum(GPP)/1000000*1800) %>%
  mutate(date = as.Date(date))

 # create daily mean and sd for EC data
flux3 <- flux2 %>%
  mutate(doy = as.numeric(strftime(date, "%j"))) %>%
  group_by(doy) %>%
  summarise(across(everything(), mean), .groups = "drop")

# load GPP from PLUMBER2 for ES-LMa 2003-2006 
plumfluxMaj<-process_nc_file("/path/to/file/ES-LMa_2004-2006_LaThuile_Flux.nc")

# Merge the two data frames (2003-2006 with 2015-2018)
maj_ec_flux_total <- bind_rows(plumfluxMaj, flux2 %>% select(date, GPP)) %>%
  filter(!(format(date, "%m") == "02" & format(date, "%d") == "29")) %>%
  mutate(date = as.numeric(strftime(date,format = "%j")))

maj_ec_flux_total <- maj_ec_flux_total[!duplicated(maj_ec_flux_total$date), ]

# calculate daily mean and sd for the EC data
maj_ec_flux_total_doymean<-maj_ec_flux_total %>%
  group_by(date) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE))))


 # Filter the QUINCY data within the date range of the EC data
ESMaGsimnormtau_meanGPP <- ESMaGsimnormtau %>%
  filter(as.Date(Date) %in% as.Date(maj_ec_flux_total$date))

 # Filter the QUINCY data within the date range of the EC data
ESMaGdev_meanGPP <- ESMaGdev %>%
  filter(as.Date(Date) %in% as.Date(maj_ec_flux_total$date)) 

 # calculate daily mean and sd for the QUINCY new version
ESMaGsimnormtau_meanGPP <- ESMaGsimnormtau_meanGPP %>%
  group_by(Time) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE))))

 # calculate daily mean and sd for the QUINCY default version
ESMaGdev_meanGPP <- ESMaGdev_meanGPP %>%
  group_by(Time) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE))))


#####################
######## GCC ########
#####################


# load gcc from phenocam
eslm1<-get_pheno_ts(site = 'eslm1', vegType = 'GR', roiID = 1000, type = '1day') %>%
  select(c('date', 'midday_gcc')) %>%
  na.omit()

range(eslm1$date)

 # load QUINCY LAI data
ESMaGsimpheno <- read.delim("path/to/quincy/dynamic/ES-MaG/output/veg_diagnostics_daily.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "sim" )

ESMaGdevpheno <- read.delim("path/to/quincy/default/ES-MaG/output/veg_diagnostics_daily.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "dev" )


eslm1_doymean <- eslm1 %>%
  select(c('date', 'midday_gcc')) %>%
  na.omit() %>%
  mutate(Date = as.Date(date),
         GCC_norm = setup::minmax(midday_gcc)) %>%
  filter(!(format(Date, "%m") == "02" & format(Date, "%d") == "29"),
         Date <= as.Date("2024-12-31")) %>% 
  select(c(3:4)) %>%
  group_by(yday(Date)) %>%
  rename(Time = 'yday(Date)') %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE))))


ESMaGsimpheno_mean <- ESMaGsimpheno %>%
  filter(Date >= as.Date("2014-06-24") & Date <= as.Date("2024-12-31")) %>%
  select(c("Time", "LAI")) %>%
  mutate(LAI_norm = setup::minmax(LAI)) %>%
  group_by(Time) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE)))
  )


ESMaGdevpheno_mean <- ESMaGdevpheno %>%
  filter(Date >= as.Date("2014-06-24") & Date <= as.Date("2024-12-31")) %>%
  select(c("Time", "LAI")) %>%
  mutate(LAI_norm = setup::minmax(LAI)) %>%
  group_by(Time) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE)))
  )



########################################################################
################################ US-Tol ################################
########################################################################

#####################
######## GPP ########
#####################

# create helper vector for assigning dates
full_date_seq <- seq(as.Date("1901-01-01"), as.Date("2025-01-01"), by = "day")
full_date_seq <- full_date_seq[!format(full_date_seq, "%m-%d") == "02-29"]

 # Load QUINCY simulation results and create Date column:
USTolsimnormtau <- read.delim("path/to/quincy/dynamic/US-TOL/output/vegfluxC_daily.txt", sep = "") %>%
  mutate(
    across("Time", round, 0),
    Date = full_date_seq[1:nrow(.)], 
    Source = "sim"
  )

USToldev <- read.delim("path/to/quincy/default/US-TOL/output/vegfluxC_daily.txt", sep = "") %>%
  mutate(
    across("Time", round, 0),
    Date = full_date_seq[1:nrow(.)], 
    Source = "sim"
  )

 # load EC flux data from FLUXNET website:

# "FLX_US-ICs_FLUXNET-CH4_HH_2014-2016_1-1.csv" 
# (from: https://fluxnet.org/data/fluxnet-ch4-community-product/)

US_ICs_EC <- read.csv("path/to/file/FLX_US-ICs_FLUXNET-CH4_HH_2014-2016_1-1.csv") %>%
  mutate(Datetime = as.POSIXct(strptime(TIMESTAMP_START, "%Y%m%d%H%M"))) %>%
  select(c("Datetime", "GPP_NT"))


US_ICs_EC_daily <- US_ICs_EC %>%
  mutate(Date = as.Date(Datetime), 
         GPP_NT = GPP_NT * 1800/1e6) %>% # convert from umol s-1 to Mol 30 min-1
  group_by(Date) %>%
  summarize(GPP_NT = sum(GPP_NT)) %>%
  na.omit()%>%  
  filter(!(format(Date, "%m") == "02" & format(Date, "%d") == "29")) 


 # match simulation with dates range of flux data 
USTolsimnormtau_IC_match <- USTolsimnormtau %>%
  filter(Date %in% US_ICs_EC_daily$Date) 

USToldev_IC_match <- USToldev %>%
  filter(Date %in% US_ICs_EC_daily$Date) 

USTolsimnormtau_IC_match_mean <- USTolsimnormtau_IC_match %>%
  group_by(Time) %>%
  summarise(GPP_mean = mean(GPP, na.rm=T),
            GPP_sd = sd(GPP, na.rm=T))

USToldev_IC_match_mean <- USToldev_IC_match %>%
  group_by(Time) %>%
  summarise(GPP_mean = mean(GPP, na.rm=T),
            GPP_sd = sd(GPP, na.rm=T))


# create helper vector
mandoy <- rep(1:365, times = 10)

US_ICs_EC_daily_mean <- US_ICs_EC_daily %>%
  mutate(GPP_NT = ifelse(GPP_NT < 0, 0, GPP_NT),
         Time = head(mandoy, nrow(.))) %>%
  group_by(Time) %>%
  summarise(GPP_mean = mean(GPP_NT, na.rm=T),
            GPP_sd = sd(GPP_NT, na.rm=T))


#####################
######## GCC ########
#####################

us_tool_gcc <- get_pheno_ts(site = 'NEON.D18.TOOL.DP1.00033', vegType = 'TN', roiID = 1000, type = '1day') %>%
  select(c('date', 'midday_gcc')) %>%
  na.omit()

USTolsimnormtau_LAI <- read.delim("path/to/quincy/dynamic/US-TOL/outpUS-TOL/outputstics_daily.txt", sep = "") %>%
  mutate(
    across("Time", round, 0),
    Date = full_date_seq[1:nrow(.)], 
    Source = "sim"
  )

USToldev_LAI <- read.delim("path/to/quincy/default/US-TOL/output/veg_diagnostics_daily.txt", sep = "") %>%
  mutate(
    across("Time", round, 0),
    Date = full_date_seq[1:nrow(.)], 
    Source = "dev"
  )

us_tool_gcc_minmax <- us_tool_gcc %>%
  filter(midday_gcc >= 0.32) %>%
  mutate(Date = as.Date(date),
         GCC_norm = setup::minmax(midday_gcc)) %>%
  filter(!(format(Date, "%m") == "02" & format(Date, "%d") == "29"))


USTolsimnormtau_LAI_minmax <- USTolsimnormtau_LAI %>%
  select(c("Time", "LAI", "Date")) %>%
  mutate(LAI_norm = setup::minmax(LAI)) %>%
  filter(Date >= as.Date("2017-01-01") & Date <= as.Date("2024-12-31")) %>%
  group_by(Time) %>%
  summarize(mean_LAI_norm = mean(LAI_norm),
            sd_LAI_norm = sd(LAI_norm)) %>%
  ungroup()

USToldev_LAI_minmax <- USToldev_LAI %>%
  select(c("Time", "LAI", "Date")) %>%
  mutate(LAI_norm = setup::minmax(LAI)) %>%
  filter(Date >= as.Date("2017-01-01") & Date <= as.Date("2024-12-31")) %>%
  group_by(Time) %>%
  summarize(mean_LAI_norm = mean(LAI_norm),
            sd_LAI_norm = sd(LAI_norm)) %>%
  ungroup()

us_tool_gcc_minmax_mean <- us_tool_gcc_minmax %>%
  group_by(yday(Date)) %>%
  summarize(mean_gcc_norm = mean(GCC_norm),
            sd_gcc_norm = sd(GCC_norm)) %>%
  ungroup() %>%
  rename(Time = 'yday(Date)')


########################################################################
################################ IE-Dri ################################
########################################################################

#####################
######## GPP ########
#####################

plumfluxIEDri<-process_nc_file("path/to/file/IE-Dri_2003-2005_LaThuile_Flux.nc")

plumfluxIEDri$doy<-strftime(plumfluxIEDri$date,format = "%j")

plumfluxIEDri$doy<-as.numeric(plumfluxIEDri$doy)
plumfluxIEDri <- plumfluxIEDri %>%
  filter(!(format(date, "%m") == "02" & format(date, "%d") == "29"))

plumfluxIEDri_doymean<-plumfluxIEDri %>%
  group_by(doy) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE))))


start_dateDri <- as.Date(min(range(plumfluxIEDri$date)))
end_date2Dri <- as.Date(max(range(plumfluxIEDri$date)))


IEDrisimnormtau <- read.delim("path/to/quincy/dynamic/IE-Dri/output/vegfluxC_daily.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "sim" )


IEDridev <- read.delim("path/to/quincy/default/IE-Dri/output/vegfluxC_daily.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "dev" )


# Filter the data within the date range
IEDrisimnormtau_meanGPP <- IEDrisimnormtau[
  IEDrisimnormtau$Date >= start_dateDri &
    IEDrisimnormtau$Date <= end_date2Dri,
]

IEDrisimnormtau_meanGPP <- IEDrisimnormtau_meanGPP %>%
  group_by(Time) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE)))
  )

IEDri_simnormtau_calcecrange <- IEDrisimnormtau[
  as.Date(IEDrisimnormtau$Date) >= start_dateDri &
    as.Date(IEDrisimnormtau$Date) <= end_date2Dri,
]

# Filter the data within the date range
IEDridev_meanGPP <- IEDridev[
  IEDridev$Date >= start_dateDri &
    IEDridev$Date <= end_date2Dri,
]

IEDridev_meanGPP <- IEDridev_meanGPP %>%
  group_by(Time) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE)))
  )

IEDri_dev_calcecrange <- IEDridev[
  IEDridev$Date >= start_dateDri &
    IEDridev$Date <= end_date2Dri,
]



########################################################################
################################# CGE ##################################
########################################################################

#####################
######## GCC ########
#####################

 # load QUINCY LAI 

CGEsimnormtau_dia <- read.delim("/path/to/quincy/dynamic/CGE_000/output/veg_diagnostics_daily.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "sim" )

CGEdev_dia <- read.delim("path/to/quincy/default/CGE_000/output/veg_diagnostics_daily.txt", sep = "") %>%
  mutate(across("Time", round, 0))%>%
  mutate(year = head(yearvec12spin, nrow(.))) %>%
  mutate(Date = paste(Time,year+1900, sep = "/"))  %>%
  mutate(Date = as.POSIXct(Date, format="%j/%Y")) %>%
  mutate( Source = "dev" )

# Assuming '1day' type for this example
CGE000 <- get_pheno_ts(site = 'gumpenstein', vegType = 'GR', roiID = 1000, type = '1day') %>%
  select(c('date', 'midday_gcc')) %>%
  na.omit()

range(CGE000$date)

CGE_doymean<-CGE000 %>%
  mutate(Date = as.Date(date),
         GCC_norm = setup::minmax(midday_gcc)) %>%
  filter(!(format(Date, "%m") == "02" & format(Date, "%d") == "29"),
         Date <= as.Date("2024-12-31")) %>% 
  select(c(3:4)) %>%
  group_by(yday(Date)) %>%
  rename(Time = 'yday(Date)') %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE))))


CGEsimnormtau_meandia <- CGEsimnormtau_dia %>%
  filter(Date >= as.Date("2018-03-24") & Date <= as.Date("2024-12-31"))

# Filter the data within the date range
CGEdev_meandia <- CGEdev_dia %>%
  filter(Date >= as.Date("2018-03-24") & Date <= as.Date("2024-12-31"))

CGEsimnormtau_meandia <- CGEsimnormtau_meandia %>%
  select(c("Time", "LAI")) %>%
  mutate(LAI_norm = setup::minmax(LAI)) %>%
  group_by(Time) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE)))
  )

CGEdev_meandia <- CGEdev_meandia %>%
  select(c("Time", "LAI")) %>%
  mutate(LAI_norm = setup::minmax(LAI)) %>%
  group_by(Time) %>%
  summarize(
    across(where(is.numeric),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE)))
  )



#################################################################################################
######################################### PLOT FIGURE 3 #########################################
#################################################################################################


ES_GPP <- ggplot() +
  # Add EC flux mean line and SD ribbon
  geom_ribbon(data = maj_ec_flux_total_doymean,
              aes(x = date, ymin = GPP_mean*12.0107 - GPP_sd*12.0107, ymax = GPP_mean*12.0107 + GPP_sd*12.0107, fill = "EC flux"),
              alpha = 0.2) +
  geom_line(data = maj_ec_flux_total_doymean,
            aes(x = date, y = GPP_mean*12.0107, color = "EC flux")) +
  
  # Add default turnover mean line and SD ribbon
  geom_ribbon(data = ESMaGdev_meanGPP,
              aes(x = Time, ymin = GPP_mean*12.0107 - GPP_sd*12.0107, ymax = GPP_mean*12.0107 + GPP_sd*12.0107, fill = "default turnover"),
              alpha = 0.2) +
  geom_line(data = ESMaGdev_meanGPP,
            aes(x = Time, y = GPP_mean*12.0107, color = "default turnover")) +
  
  # Add new turnover mean line and SD ribbon
  geom_ribbon(data = ESMaGsimnormtau_meanGPP,
              aes(x = Time, ymin = GPP_mean*12.0107 - GPP_sd*12.0107, ymax = GPP_mean*12.0107 + GPP_sd*12.0107, fill = "dynamic turnover"),
              alpha = 0.2) +
  geom_line(data = ESMaGsimnormtau_meanGPP,
            aes(x = Time, y = GPP_mean*12.0107, color = "dynamic turnover")) +
  
  labs(
    title = "ES-LMa mean GPP for 2003-2006 and 2015-2018",
    x = "Day of Year",
    y = "GPP (g C m<sup>-2</sup> day<sup>-1</sup>)")+
  
  # Define color and fill scales
  scale_color_manual(
    name = "Data",
    breaks = c("default turnover", "dynamic turnover", "EC flux"),
    values = c("default turnover" = cbbPalette[6], "dynamic turnover" = "#ff641e", "EC flux" = cbbPalette[4])
  ) +
  scale_fill_manual(
    name = "Data",
    breaks = c("default turnover", "dynamic turnover", "EC flux"),
    values = c("default turnover" = cbbPalette[6], "dynamic turnover" = "#ff641e", "EC flux" = cbbPalette[4])
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.y = element_markdown(size = 36, face = 'plain'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) +
  scale_y_continuous(limits = c(0, 9.5), breaks = seq(0, 9.5, 2))+
  scale_x_continuous(limits = c(0, 366), breaks = seq(0, 365, 50)); print(ES_GPP)




ESMaG_gcc <- ggplot() +   # Add SD ribbon for PhenoCam
  geom_ribbon(data = eslm1_doymean,
              aes(x = Time,
                  ymin = GCC_norm_mean - GCC_norm_sd,
                  ymax = GCC_norm_mean + GCC_norm_sd,
                  fill = "PhenoCam GCC"), alpha = 0.2) +  # Adjust alpha for transparency
  
  # Add EC flux mean line and SD ribbon for PhenoCam
  geom_line(data = eslm1_doymean,
            aes(x = Time, y = GCC_norm_mean, color = "PhenoCam GCC")) +
  
  # Add SD ribbon for PhenoCam
  geom_ribbon(data = ESMaGdevpheno_mean,
              aes(x = Time,
                  ymin = LAI_norm_mean - LAI_norm_sd,
                  ymax = LAI_norm_mean + LAI_norm_sd,
                  fill = "LAI default turnover"), alpha = 0.2) +  # Adjust alpha for transparency
  
  # Add mean line for default turnover
  geom_line(data = ESMaGdevpheno_mean,
            aes(x = Time, y = LAI_norm_mean, color = "LAI default turnover")) +
  
  # Add SD ribbon for default turnover (if needed)
  geom_ribbon(data = ESMaGsimpheno_mean,
              aes(x = Time,
                  ymin = LAI_norm_mean - LAI_norm_sd,
                  ymax = LAI_norm_mean + LAI_norm_sd,
                  fill = "LAI dynamic turnover"), alpha = 0.2) +  # Adjust alpha for transparency
  
  # Add mean line for dynamic turnover
  geom_line(data = ESMaGsimpheno_mean,
            aes(x = Time, y = LAI_norm_mean, color = "LAI dynamic turnover")) +
  labs(
    title = "ES-LMa mean LAI and GCC for 2014-2024",
    x = "Day of Year",
    y = "")+
  
  # Define color and fill scales
  scale_color_manual(
    name = "Data",
    breaks = c("LAI default turnover", "LAI dynamic turnover", "PhenoCam GCC"),
    values = c("LAI default turnover" = cbbPalette[6], "LAI dynamic turnover" = "#ff641e", "PhenoCam GCC" = cbbPalette[4])
  ) +
  scale_fill_manual(
    name = "Data",
    breaks = c("LAI default turnover", "LAI dynamic turnover", "PhenoCam GCC"),
    values = c("LAI default turnover" = cbbPalette[6], "LAI dynamic turnover" = "#ff641e", "PhenoCam GCC" = cbbPalette[4])
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.y = element_markdown(size = 36, face = 'plain'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) +
  coord_cartesian(ylim = c(0, .85))+
  scale_y_continuous(limits = c(0, .9), breaks = seq(0, .9, 0.1)) +
  scale_x_continuous(limits = c(0, 365), breaks = seq(0, 365, 50)) ; print(ESMaG_gcc)



##########################################################################################
#######################################  TOOLIK ##########################################
##########################################################################################


US_GPP <- ggplot()+
  geom_ribbon(data=US_ICs_EC_daily_mean,aes(x = Time,
                                            ymin = (GPP_mean*12.0107 - GPP_sd*12.0107),
                                            ymax = (GPP_mean*12.0107 + GPP_sd*12.0107), fill = "EC flux"), alpha=.3)+
  geom_line(data=US_ICs_EC_daily_mean, aes(x=Time, y = GPP_mean*12.0107, color = "EC flux"))+
  
  geom_ribbon(data=USTolsimnormtau_IC_match_mean,aes(x = Time,
                                                     ymin = (GPP_mean*12.0107 - GPP_sd*12.0107),
                                                     ymax = (GPP_mean*12.0107 + GPP_sd*12.0107), fill = "dynamic turnover"), alpha=.3)+
  
  geom_line(data=USTolsimnormtau_IC_match_mean,aes(x=Time, y = (GPP_mean*12.0107), color = "dynamic turnover"))+
  
  geom_ribbon(data=USToldev_IC_match_mean,aes(x = Time,
                                              ymin = (GPP_mean*12.0107 - GPP_sd*12.0107),
                                              ymax = (GPP_mean*12.0107 + GPP_sd*12.0107), fill = "default turnover"), alpha=.3)+
  
  geom_line(data=USToldev_IC_match_mean,aes(x=Time, y = (GPP_mean*12.0107), color = "default turnover"))+
  labs(
    title = "US-Tol mean GPP for 2014-2016",
    x = "Day of Year",
    y = "GPP (g C m<sup>-2</sup> day<sup>-1</sup>)")+
  
  # Define color and fill scales
  scale_color_manual(
    name = "Data",
    breaks = c("default turnover", "dynamic turnover", "EC flux"),
    values = c("default turnover" = cbbPalette[6], "dynamic turnover" = "#ff641e", "EC flux" = cbbPalette[4])
  ) +
  scale_fill_manual(
    name = "Data",
    breaks = c("default turnover", "dynamic turnover", "EC flux"),
    values = c("default turnover" = cbbPalette[6], "dynamic turnover" = "#ff641e", "EC flux" = cbbPalette[4])
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.y = element_markdown(size = 36, face = 'plain'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) +
  coord_cartesian(ylim = c(0, 5))+
  scale_y_continuous(limits = c(-1, 5), breaks = seq(-1, 5, 1)) +
  scale_x_continuous(limits = c(0, 365), breaks = seq(0, 365, 50)); print(US_GPP)


USTol_gcc <- ggplot() +
  geom_line(data = us_tool_gcc_minmax_mean, aes(x = Time, y = mean_gcc_norm, color = "GCC")) +
  geom_ribbon(data = us_tool_gcc_minmax_mean, aes(x = Time, ymin = mean_gcc_norm - sd_gcc_norm, ymax = mean_gcc_norm + sd_gcc_norm), 
              fill = cbbPalette[4], alpha = 0.3) +
  
  geom_line(data = USToldev_LAI_minmax, aes(x = Time, y = mean_LAI_norm, color = "LAI dynamic turnover" )) +
  geom_ribbon(data = USToldev_LAI_minmax, aes(x = Time, ymin = mean_LAI_norm - sd_LAI_norm, ymax = mean_LAI_norm + sd_LAI_norm), 
              fill = "#ff641e", alpha = 0.3) +
  
  geom_line(data = USTolsimnormtau_LAI_minmax, aes(x = Time, y = mean_LAI_norm, color = "LAI default turnover"  )) +
  geom_ribbon(data = USTolsimnormtau_LAI_minmax, aes(x = Time, ymin = mean_LAI_norm - sd_LAI_norm, ymax = mean_LAI_norm + sd_LAI_norm), 
              fill = cbbPalette[6], alpha = 0.3) +
  labs(
    title = "US-Tol mean LAI and GCC for 2017-2024",
    x = "Day of Year",
    y = "")+
  
  # Define color and fill scales
  scale_color_manual(
    name = "Data",
    breaks = c("LAI default turnover", "LAI dynamic turnover", "GCC"),
    values = c("LAI default turnover" = cbbPalette[6], "LAI dynamic turnover" = "#ff641e", "GCC" = cbbPalette[4])
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.y = element_markdown(size = 36, face = 'plain'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) +
  coord_cartesian(ylim = c(0, .95))+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(limits = c(0, 365), breaks = seq(0, 365, 50)) ; print(USTol_gcc)


##########################################################################################
####################################  CLIMGRASS ##########################################
##########################################################################################


CGE_gcc <- ggplot() +   # Add SD ribbon for PhenoCam
  geom_ribbon(data = CGE_doymean,
              aes(x = Time,
                  ymin = GCC_norm_mean - GCC_norm_sd,
                  ymax = GCC_norm_mean + GCC_norm_sd,
                  fill = "PhenoCam GCC"), alpha = 0.2) +  # Adjust alpha for transparency
  
  # Add EC flux mean line and SD ribbon for PhenoCam
  geom_line(data = CGE_doymean,
            aes(x = Time, y = GCC_norm_mean, color = "PhenoCam GCC")) +
  
  # Add SD ribbon for PhenoCam
  geom_ribbon(data = CGEdev_meandia,
              aes(x = Time,
                  ymin = LAI_norm_mean - LAI_norm_sd,
                  ymax = LAI_norm_mean + LAI_norm_sd,
                  fill = "LAI default turnover"), alpha = 0.2) +  # Adjust alpha for transparency
  
  # Add mean line for default turnover
  geom_line(data = CGEdev_meandia,
            aes(x = Time, y = LAI_norm_mean, color = "LAI default turnover")) +
  
  # Add SD ribbon for default turnover (if needed)
  geom_ribbon(data = CGEsimnormtau_meandia,
              aes(x = Time,
                  ymin = LAI_norm_mean - LAI_norm_sd,
                  ymax = LAI_norm_mean + LAI_norm_sd,
                  fill = "LAI dynamic turnover"), alpha = 0.2) +  # Adjust alpha for transparency
  
  # Add mean line for dynamic turnover
  geom_line(data = CGEsimnormtau_meandia,
            aes(x = Time, y = LAI_norm_mean, color = "LAI dynamic turnover")) +
  labs(
    title = "CGE mean LAI and GCC for 2018-2024",
    x = "Day of Year",
    y = "")+
  
  # Define color and fill scales
  scale_color_manual(
    name = "Data",
    breaks = c("LAI default turnover", "LAI dynamic turnover", "PhenoCam GCC"),
    values = c("LAI default turnover" = cbbPalette[6], "LAI dynamic turnover" = "#ff641e", "PhenoCam GCC" = cbbPalette[4])
  ) +
  scale_fill_manual(
    name = "Data",
    breaks = c("LAI default turnover", "LAI dynamic turnover", "PhenoCam GCC"),
    values = c("LAI default turnover" = cbbPalette[6], "LAI dynamic turnover" = "#ff641e", "PhenoCam GCC" = cbbPalette[4])
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.y = element_markdown(size = 36, face = 'plain'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) + 
  coord_cartesian(ylim = c(0, .95)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(limits = c(0, 365), breaks = seq(0, 365, 50)) ; print(CGE_gcc)



##########################################################################################
######################################  DRIPSEY ##########################################
##########################################################################################


IE_GPP <- ggplot() +
  # Add EC flux mean line and SD ribbon
  geom_ribbon(data = plumfluxIEDri_doymean,
              aes(x = doy, ymin = GPP_mean*12.0107 - GPP_sd*12.0107, ymax = GPP_mean*12.0107 + GPP_sd*12.0107, fill = "EC flux"),
              alpha = 0.2) +
  geom_line(data = plumfluxIEDri_doymean,
            aes(x = doy, y = GPP_mean*12.0107, color = "EC flux")) +
  
  # Add default turnover mean line and SD ribbon
  geom_ribbon(data = IEDridev_meanGPP,
              aes(x = Time, ymin = GPP_mean*12.0107 - GPP_sd*12.0107, ymax = GPP_mean*12.0107 + GPP_sd*12.0107, fill = "default turnover"),
              alpha = 0.2) +
  geom_line(data = IEDridev_meanGPP,
            aes(x = Time, y = GPP_mean*12.0107, color = "default turnover")) +
  
  # Add dynamic turnover mean line and SD ribbon
  geom_ribbon(data = IEDrisimnormtau_meanGPP,
              aes(x = Time, ymin = GPP_mean*12.0107 - GPP_sd*12.0107, ymax = GPP_mean*12.0107 + GPP_sd*12.0107, fill = "dynamic turnover"),
              alpha = 0.2) +
  geom_line(data = IEDrisimnormtau_meanGPP,
            aes(x = Time, y = GPP_mean*12.0107, color = "dynamic turnover")) +
  
  labs(
    title = "IE-Dri mean GPP for 2003-2005",
    x = "Day of Year",
    y = "GPP (g C m<sup>-2</sup> day<sup>-1</sup>)")+
  
  # Define color and fill scales
  scale_color_manual(
    name = "Data",
    breaks = c("default turnover", "dynamic turnover", "EC flux"),
    values = c("default turnover" = cbbPalette[6], "dynamic turnover" = "#ff641e", "EC flux" = cbbPalette[4])
  ) +
  scale_fill_manual(
    name = "Data",
    breaks = c("default turnover", "dynamic turnover", "EC flux"),
    values = c("default turnover" = cbbPalette[6], "dynamic turnover" = "#ff641e", "EC flux" = cbbPalette[4])
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_markdown(size = 26, face = 'bold'),
    axis.title.y = element_markdown(size = 36, face = 'plain'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   # x-axis tick labels
    axis.text.y = element_text(size = 36, face = 'plain'),   # y-axis tick labels
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) +
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 15, 2))+
  scale_x_continuous(limits = c(0, 366), breaks = seq(0, 365, 50)); print(IE_GPP)


##########################################################################################
######################################## CREATE ##########################################
##########################################################################################

legendIEDri <- get_legend(IE_GPP)

ES_GPP2 <- ES_GPP + 
  annotate("text", x = -Inf, y = Inf, label = "a)", hjust = -0.2, vjust = 1.2, size = 12, fontface = "bold") + 
  theme(axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(size = 1),
        legend.position = "none")

US_GPP2 <- US_GPP + 
  annotate("text", x = -Inf, y = Inf, label = "c)", hjust = -0.2, vjust = 1.2, size = 12, fontface = "bold") + 
  theme(axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(size = 1),
        legend.position = "none")

IE_GPP2 <- IE_GPP + 
  annotate("text", x = -Inf, y = Inf, label = "e)", hjust = -0.2, vjust = 1.2, size = 12, fontface = "bold") + 
  theme(axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(size = 1),
        legend.position = "none")

bigfour <- ggarrange(ES_GPP2, 
                     US_GPP2, 
                     IE_GPP2,
                     legendIEDri, nrow = 4, heights = c(1,1,1,.2)) ; print(bigfour)

legendUSTOL_GCC <- get_legend(USTol_gcc)

ESMaG_gcc2 <- ESMaG_gcc + 
  annotate("text", x = -Inf, y = Inf, label = "b)", hjust = -0.2, vjust = 1.2, size = 12, fontface = "bold") + 
  theme(axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(size = 1),
        legend.position = "none")

USTol_gcc2 <- USTol_gcc + 
  annotate("text", x = -Inf, y = Inf, label = "d)", hjust = -0.2, vjust = 1.2, size = 12, fontface = "bold") + 
  theme(axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(size = 1),
        legend.position = "none")

CGE_gcc2 <- CGE_gcc + 
  annotate("text", x = -Inf, y = Inf, label = "f)", hjust = -0.2, vjust = 1.2, size = 12, fontface = "bold") + 
  theme(axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(size = 1),
        legend.position = "none")

bigfour_GCC <- ggarrange(ESMaG_gcc2, 
                         USTol_gcc2, 
                         CGE_gcc2, 
                         legendUSTOL_GCC,  nrow = 4, heights = c(1,1,1,.2)) ; print(bigfour_GCC)

bigfour_all_GPP_GCC <- ggarrange(bigfour, bigfour_GCC, ncol = 2, widths = c(1,1)); print(bigfour_all_GPP_GCC)  # Adjust heights


ggsave(
  "path/to/figure.png",
  plot = bigfour_all_GPP_GCC,
  scale = 6,
  width = 1000,
  height = 1600,
  units = "px")




#################################################################################################
######################################### STATS TABLE 4 #########################################
#################################################################################################



r_squared <- function(observed, predicted) {
  
  model = lm(observed~predicted) 
  summary(model)
  summary(model)$r.squared 
}


sites_pap <- c("cge_dev", "cge_new", "esmag_dev", "esmag_new", "ustol_dev", "ustol_new", "iedri_dev", "iedri_new")

df <- data.frame(
  mae = rep(NA, 8),
  rmse = rep(NA, 8),
  nrmse = rep(NA, 8),
  Pearson = rep(NA, 8),
  r2 = rep(NA, 8),
  row.names = sites_pap
); print(df)

####################################################################################################
############################################## IE-DRI ##############################################
####################################################################################################

IEobs <- plumfluxIEDri_doymean$GPP_mean[-60]
IEdev <- IEDridev_meanGPP$GPP_mean
IEnew <- IEDrisimnormtau_meanGPP$GPP_mean

######### RMSE AND NRMSE #####################
# dev
mae_dev_IE <- hydroGOF::mae(obs = IEobs,sim =  IEdev, na.rm=T)
print(mae_dev_IE) 
df["iedri_dev", "mae"] <- mae_dev_IE

# dev
rmse_dev_IE <- hydroGOF::rmse(obs = IEobs,sim =  IEdev, na.rm=T)
#rmse_dev_IE_uncertainty <- rmse_se(IEobs, IEdev)
print(rmse_dev_IE) 
df["iedri_dev", "rmse"] <- rmse_dev_IE
#print(rmse_dev_IE_uncertainty)

nrmse_dev_IE <- INDperform::nrmse(obs = IEobs, pred = IEdev)
#nrmse_dev_IE_uncertainty <-nrmse_se(IEobs, IEdev, method = "range" )
print(nrmse_dev_IE)
#print(nrmse_dev_IE_uncertainty)
df["iedri_dev", "nrmse"] <- nrmse_dev_IE

# new
mae_new_IE <- hydroGOF::mae(obs = IEobs,sim =  IEnew, na.rm=T)
print(mae_new_IE) 
df["iedri_new", "mae"] <- mae_new_IE

rmse_new_IE <- hydroGOF::rmse(obs = IEobs,sim =  IEnew, na.rm=T)
#rmse_new_IE_uncertainty <- rmse_se(IEobs, IEnew)
print(rmse_new_IE)
#print(rmse_new_IE_uncertainty)
df["iedri_new", "rmse"] <- rmse_new_IE


nrmse_new_IE <- INDperform::nrmse(obs = IEobs, pred = IEnew)
#nrmse_new_IE_uncertainty <-nrmse_se(IEobs, IEnew, method = "range" )
print(nrmse_new_IE)
#print(nrmse_new_IE_uncertainty)
df["iedri_new", "nrmse"] <- nrmse_new_IE


################## pearson's ################## 

# dev
iedevcor <- cor(IEdev, IEobs, method = "pearson", use = "complete.obs")
df["iedri_dev", "Pearson"] <- iedevcor

# new
ienewcor <- cor(IEnew, IEobs, method = "pearson", use = "complete.obs")
df["iedri_new", "Pearson"] <- ienewcor


################## R^2 ##################  

# dev
iedevr2 <- r_squared(observed = plumfluxIEDri$GPP,predicted = IEDri_dev_calcecrange$GPP )
df["iedri_dev", "r2"] <- iedevr2


# new
ienewr2 <- r_squared(observed = plumfluxIEDri$GPP, predicted = IEDri_simnormtau_calcecrange$GPP)
df["iedri_new", "r2"] <- ienewr2


####################################################################################################
############################################## US-Tol ##############################################
####################################################################################################

USobs <- US_ICs_EC_daily_mean$GPP_mean
USdev <- USToldev_IC_match_mean$GPP_mean
USnew <- USTolsimnormtau_IC_match_mean$GPP_mean

######### RMSE AND NRMSE #####################
# dev
mae_dev_US <- hydroGOF::mae(obs = USobs,sim = USdev, na.rm=T)
print(mae_dev_US) 
df["ustol_dev", "mae"] <- mae_dev_US

rmse_dev_US <-hydroGOF::rmse(obs = USobs,sim =  USdev, na.rm=T)
#rmse_dev_US_uncertainty <- rmse_se(USobs, USdev)
print(rmse_dev_US)
#print(rmse_dev_US_uncertainty)
df["ustol_dev", "rmse"] <- rmse_dev_US

nrmse_dev_US <- INDperform::nrmse(obs = USobs, pred = USdev)
#nrmse_dev_US_uncertainty <-nrmse_se(USobs, USdev, method = "range" )
print(nrmse_dev_US)
#print(nrmse_dev_US_uncertainty)
df["ustol_dev", "nrmse"] <- nrmse_dev_US

# new
mae_new_US <- hydroGOF::mae(obs = USobs,sim =  USnew, na.rm=T)
print(mae_new_US) 
df["ustol_new", "mae"] <- mae_new_US

rmse_new_US <- hydroGOF::rmse(obs = USobs,sim =  USnew, na.rm=T)
#rmse_new_US_uncertainty <- rmse_se(USobs, USnew)
print(rmse_new_US)
#print(rmse_new_US_uncertainty)
df["ustol_new", "rmse"] <- rmse_new_US

nrmse_new_US <- INDperform::nrmse(obs = USobs, pred = USnew)
#nrmse_new_US_uncertainty <-nrmse_se(USobs, USnew, method = "range" )
print(nrmse_new_US)
#print(nrmse_new_US_uncertainty)
df["ustol_new", "nrmse"] <- nrmse_new_US

################## pearson's ################## 

# dev
usdevcor <- cor(USdev, USobs, method = "pearson", use = "complete.obs")
df["ustol_dev", "Pearson"] <- usdevcor
# new
usnewcor <- cor(USnew, USobs, method = "pearson", use = "complete.obs")
df["ustol_new", "Pearson"] <- usnewcor

################## R^2 ##################  

# dev
usdevr2 <- r_squared(predicted = USdev, observed = USobs)
df["ustol_dev", "r2"] <- usdevr2

# new
usnewr2 <- r_squared(predicted = USnew, observed = USobs)
df["ustol_new", "r2"] <- usnewr2



####################################################################################################
############################################## ES-MaG ##############################################
####################################################################################################

ESobs <- maj_ec_flux_total_doymean$GPP_mean[-60]
ESdev <- ESMaGdev_meanGPP$GPP_mean
ESnew <- ESMaGsimnormtau_meanGPP$GPP_mean

######### RMSE AND NRMSE #####################
# dev
mae_dev_ES <- hydroGOF::mae(obs = ESobs,sim =  ESdev, na.rm=T)
print(mae_dev_ES) 
df["esmag_dev", "mae"] <- mae_dev_ES

rmse_dev_ES <-hydroGOF::rmse(obs = ESobs,sim =  ESdev, na.rm=T)
#rmse_dev_ES_uncertainty <- rmse_se(ESobs, ESdev)
print(rmse_dev_ES)
#print(rmse_dev_ES_uncertainty)
df["esmag_dev", "rmse"] <- rmse_dev_ES


nrmse_dev_ES <- INDperform::nrmse(obs = ESobs, pred = ESdev)
#nrmse_dev_ES_uncertainty <-nrmse_se(ESobs, ESdev, method = "range" )
print(nrmse_dev_ES)
#print(nrmse_dev_ES_uncertainty)
df["esmag_dev", "nrmse"] <- nrmse_dev_ES

# new
mae_new_ES <- hydroGOF::mae(obs = ESobs,sim =  ESnew, na.rm=T)
print(mae_new_ES) 
df["esmag_new", "mae"] <- mae_new_ES

rmse_new_ES <- hydroGOF::rmse(obs = ESobs,sim =  ESnew, na.rm=T)
#rmse_new_ES_uncertainty <- rmse_se(ESobs, ESnew)
print(rmse_new_ES)
#print(rmse_new_ES_uncertainty)
df["esmag_new", "rmse"] <- rmse_new_ES


nrmse_new_ES <- INDperform::nrmse(obs = ESobs, pred = ESnew)
#nrmse_new_ES_uncertainty <-nrmse_se(ESobs, ESnew, method = "range" )
print(nrmse_new_ES)
#print(nrmse_new_ES_uncertainty)
df["esmag_new", "nrmse"] <- nrmse_new_ES


################## pearson's ################## 

# dev
ESdevcor <- cor(ESdev, ESobs, method = "pearson", use = "complete.obs")
df["esmag_dev", "Pearson"] <- ESdevcor

# new
ESnewcor <- cor(ESnew, ESobs, method = "pearson", use = "complete.obs")
df["esmag_new", "Pearson"] <- ESnewcor


################## R^2 ##################  

# dev
ESdevr2 <- r_squared(observed = ESobs, predicted = ESdev)
df["esmag_dev", "r2"] <- ESdevr2

# new
ESnewr2 <- r_squared(predicted = ESnew, observed = ESobs)
df["esmag_new", "r2"] <- ESnewr2


####################################################################################################
############################################## CGE000 ##############################################
####################################################################################################

CGobs <- (CGE_doymean$GCC_norm_mean[1:365])
CGdev <- (CGEdev_meandia$LAI_norm_mean)
CGnew <- (CGEsimnormtau_meandia$LAI_norm_mean)

######### RMSE AND NRMSE #####################
# dev
mae_dev_CG <- hydroGOF::mae(obs = CGobs,sim =  CGdev, na.rm=T)
print(mae_dev_CG) 
df["cge_dev", "mae"] <- mae_dev_CG

rmse_dev_CG <- hydroGOF::rmse(obs = CGobs,sim =  CGdev, na.rm=T)
#rmse_dev_CG_uncertainty <- rmse_se(CGobs, CGdev)
print(rmse_dev_CG)
#print(rmse_dev_CG_uncertainty)
df["cge_dev", "rmse"] <- rmse_dev_CG


nrmse_dev_CG <- INDperform::nrmse(obs = CGobs, pred = CGdev)
#nrmse_dev_CG_uncertainty <-nrmse_se(CGobs, CGdev, method = "range" )
print(nrmse_dev_CG)
#print(nrmse_dev_CG_uncertainty)
df["cge_dev", "nrmse"] <- nrmse_dev_CG

# new
mae_new_CG <- hydroGOF::mae(obs = CGobs,sim =  CGnew, na.rm=T)
print(mae_new_CG) 
df["cge_new", "mae"] <- mae_new_CG

rmse_new_CG <- hydroGOF::rmse(obs = CGobs,sim =  CGnew, na.rm=T)
#rmse_new_CG_uncertainty <- rmse_se(CGobs, CGnew)
print(rmse_new_CG)
#print(rmse_new_CG_uncertainty)
df["cge_new", "rmse"] <- rmse_new_CG


nrmse_new_CG <- INDperform::nrmse(obs = CGobs, pred = CGnew)
#nrmse_new_CG_uncertainty <-nrmse_se(CGobs, CGnew, method = "range" )
print(nrmse_new_CG)
#print(nrmse_new_CG_uncertainty)
df["cge_new", "nrmse"] <- nrmse_new_CG


################## pearson's ################## 

# dev
CGdevcor <- cor(CGdev, CGobs, method = "pearson", use = "complete.obs")
df["cge_dev", "Pearson"] <- CGdevcor

# new
CGnewcor <- cor(CGnew, CGobs, method = "pearson", use = "complete.obs")
df["cge_new", "Pearson"] <- CGnewcor

################## R^2 ##################  

# dev
CGdevr2 <- r_squared(predicted = CGdev, observed = CGobs)
df["cge_dev", "r2"] <- CGdevr2

# new
CGnewr2 <- r_squared(predicted = CGnew, observed = CGobs)
df["cge_new", "r2"] <- CGnewr2

print(df)



#### GCC stats ES-LMa
hydroGOF::mae(obs = eslm1_doymean$GCC_norm_mean[1:365], sim = ESMaGsimpheno_mean$LAI_norm_mean)
hydroGOF::rmse(obs = eslm1_doymean$GCC_norm_mean[1:365], sim = ESMaGsimpheno_mean$LAI_norm_mean)
hydroGOF::nrmse(obs = eslm1_doymean$GCC_norm_mean[1:365], sim = ESMaGsimpheno_mean$LAI_norm_mean, norm = 'maxmin')

cor(eslm1_doymean$GCC_norm_mean[1:365], ESMaGsimpheno_mean$LAI_norm_mean)
summary(lm(eslm1_doymean$GCC_norm_mean[1:365]~ESMaGsimpheno_mean$LAI_norm_mean))$adj.r.squared 

hydroGOF::mae(obs = eslm1_doymean$GCC_norm_mean[1:365], sim = ESMaGdevpheno_mean$LAI_norm_mean)
hydroGOF::rmse(obs = eslm1_doymean$GCC_norm_mean[1:365], sim = ESMaGdevpheno_mean$LAI_norm_mean)
hydroGOF::nrmse(obs = eslm1_doymean$GCC_norm_mean[1:365], sim = ESMaGdevpheno_mean$LAI_norm_mean, norm = 'maxmin')

cor(eslm1_doymean$GCC_norm_mean[1:365], ESMaGdevpheno_mean$LAI_norm_mean)
summary(lm(eslm1_doymean$GCC_norm_mean[1:365]~ESMaGdevpheno_mean$LAI_norm_mean))$adj.r.squared 


#### GCC stats US-Tol
hydroGOF::mae(obs = us_tool_gcc_minmax_mean$mean_gcc_norm[1:365], sim = USTolsimnormtau_LAI_minmax$mean_LAI_norm)
hydroGOF::rmse(obs = us_tool_gcc_minmax_mean$mean_gcc_norm[1:365], sim = USTolsimnormtau_LAI_minmax$mean_LAI_norm)
hydroGOF::nrmse(obs = us_tool_gcc_minmax_mean$mean_gcc_norm[1:365], sim = USTolsimnormtau_LAI_minmax$mean_LAI_norm, norm = 'maxmin')

#cor(us_tool_gcc_minmax_mean$mean_gcc_norm[1:365], USTolsimnormtau_LAI_minmax$mean_LAI_norm, use = "complete.obs")
summary(lm(us_tool_gcc_minmax_mean$mean_gcc_norm[1:365] ~ USTolsimnormtau_LAI_minmax$mean_LAI_norm))$adj.r.squared 

hydroGOF::mae(obs = us_tool_gcc_minmax_mean$mean_gcc_norm[1:365], sim = USToldev_LAI_minmax$mean_LAI_norm)
hydroGOF::rmse(obs = us_tool_gcc_minmax_mean$mean_gcc_norm[1:365], sim = USToldev_LAI_minmax$mean_LAI_norm)
hydroGOF::nrmse(obs = us_tool_gcc_minmax_mean$mean_gcc_norm[1:365], sim = USToldev_LAI_minmax$mean_LAI_norm, norm = 'maxmin')

#cor(us_tool_gcc_minmax_mean$mean_gcc_norm[1:365], USToldev_LAI_minmax$mean_LAI_norm, use = "complete.obs")
summary(lm(us_tool_gcc_minmax_mean$mean_gcc_norm[1:365] ~ USToldev_LAI_minmax$mean_LAI_norm))$adj.r.squared 





