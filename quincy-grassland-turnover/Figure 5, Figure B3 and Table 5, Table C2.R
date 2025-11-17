

# load libraries
library(tidyverse)
library(terra)
library(ggtext)
library(phenofit)
library(fs)

# create color blind friendly palette and define shapes for figure
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
shapes <- c(0:6, 15:20, 21:25)


########################################################################
################################ SETUP #################################
########################################################################

# process all NC files and aggregate GPP to daily

process_nc_data <- function(nc_obj) {
  
  file_name <- basename(nc_obj$filename)
  site_name <- strsplit(file_name, "_")[[1]][1]  
  year <- substr(strsplit(file_name, "_")[[1]][2], 1, 4)
  
  # Replace dashes with underscores
  site_name <- gsub("-", "_", site_name)
  
  # Get variables
  GPP <- ncvar_get(nc_obj, "GPP")
  time <- ncvar_get(nc_obj, "time")
  
  # time conversion
  origin_date <- as.POSIXct(paste0(year, "-01-01 00:00:00"), tz = "CET")
  date <- origin_date + time
  
  data <- data.frame(date = date, GPP = GPP)
  data$GPP_umol_s <- (data$GPP * 1e6) / 86400
  
  # Aggregate half-hourly GPP to daily GPP (mol/day)
  data$date <- as.Date(data$date)
  data <- aggregate(GPP ~ date, data = data, FUN = function(x) sum(x * 1800) / 1e6)
  
  return(list(site = site_name, data = data))
}


########################################################################################################
################################ LOAD GPP FOR ALL PLUMBER2 GRASS SITES #################################
########################################################################################################

# PLUMBER2 scripts for obtaining flux and met data: https://github.com/aukkola/PLUMBER2

plumber2sites <- read.delim(file = "/path/to/plumber2/siteinfo/all_sites_list_t.dat", sep = "")

# List of Site IDs
site_ids <- c("ES-MaG", "CGE_000", "AT-Neu", "AU-Cpr", "AU-DaP", "AU-DaS", "AU-Dry",
              "AU-Emr", "AU-GWW", "AU-Otw", "AU-Rig", "AU-Sam", "AU-Stp", "AU-TTE",
              "AU-Ync", "BW-Ma1", "CA-NS6", "CA-NS7", "CA-SF3", "CH-Cha", "CH-Fru",
              "CH-Oe1", "CN-Cng", "CN-Dan", "CN-Du2", "CN-HaM", "CZ-wet", "DE-Gri",
              "DE-SfN", "DK-ZaH", "ES-LgS", "ES-VDA", "FI-Kaa", "FI-Lom", "FR-Lq2",
              "IE-Dri", "IT-Amp", "IT-MBo", "NL-Ca1", "NL-Hor", "PT-Mi2", "RU-Che",
              "SD-Dem", "SE-Deg", "US-AR1", "US-AR2", "US-Aud", "US-Cop", "US-FPe",
              "US-Goo", "US-Los", "US-Myb", "US-SRG", "US-Tw4", "US-Var", "US-Whs",
              "US-Wkg", "ZA-Kru")


site_prefixes <- substr(site_ids, 1, 6)
nc_files <- list.files(path = "/path/to/nc/files/", pattern = "\\.nc$", full.names = TRUE)

matching_files <- nc_files[sapply(nc_files, function(x) substr(basename(x), 1, 6) %in% site_prefixes)]
nc_data <- lapply(matching_files, nc_open)

processed_data <- lapply(nc_data, process_nc_data)

names(processed_data) <- sapply(processed_data, function(x) x$site)

# Combine all processed data into a single data frame
combined_data <- do.call(rbind, lapply(processed_data, function(x) {
  x$data$site <- x$site  # Add site column
  return(x$data)
}))

######################################################################################################
################################### LOAD DYNAMIC QUINCY GPP DATA #####################################
######################################################################################################

setwd("/path/to/quincy/dynamic/")

site_folders <- dir(pattern = "^[A-Z]{2}-[A-Za-z0-9]+")

combined_datasim2 <- list()

# read all vegflux C files (contains GPP data in QUINCY):
for (site in site_folders) {
  file_path <- file.path(site, "output", "vegfluxC_daily.txt")
  
  if (file.exists(file_path)) {
    site_data <- read.delim(file_path, sep = "")
    
    site_data$Site <- site
    
    combined_datasim2[[site]] <- site_data
  } else {
    warning(paste("File not found for site:", site))
  }
}

# Combine all data frames into one big data frame
final_data2 <- do.call(rbind, combined_datasim2)

# helper vector for years
yearvec_template <- rep(1:124, each = 365)
yearvec_all_sites <- rep(yearvec_template, times = 57)

final_data_cut <- final_data2 %>%
  mutate(year = yearvec_all_sites)

# Add Date to final_data2, resetting for each site
final_data_cut <- final_data_cut %>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(across("Time", round, 0))%>%
  mutate( # Adjust length to match site data
    Date = paste(Time, year + 1900, sep = "/"),  
    Date = as.POSIXct(Date, format = "%j/%Y")    
  )

final_data_cut <- final_data_cut %>%
  mutate(Site = gsub("-", "_", Site))  

# combine EC data with dynamic turnover data
matched_data <- final_data_cut %>%
  inner_join(combined_data, by = c("Site" = "site", "Date" = "date"))

complete_years_data_sim_gpp <- matched_data %>%
  group_by(Site, year) %>%
  summarize(days_in_year = n(), .groups = "drop") %>%  
  filter(days_in_year == 365 | days_in_year == 366) %>%
  inner_join(matched_data, by = c("Site", "year")) %>% 
  select(Site, Time, Date, year, GPP.x, GPP.y) %>%
  mutate(GPP.y = ifelse(GPP.y < 0, 0, GPP.y))

# calculate annual GPP sum
annualsumGPP_sim <- complete_years_data_sim_gpp %>%
  group_by(year, Site) %>%
  summarize(
    total_GPP.x = sum(GPP.x, na.rm = TRUE),
    total_GPP.y = sum(GPP.y, na.rm = TRUE),
    .groups = "drop" 
  ) %>%
  group_by(Site) %>%
  summarize(
    mean_annual_GPP.x = mean(total_GPP.x, na.rm = TRUE),
    sd_annual_GPP.x = sd(total_GPP.x, na.rm = TRUE), 
    mean_annual_GPP.y = mean(total_GPP.y, na.rm = TRUE),
    sd_annual_GPP.y = sd(total_GPP.y, na.rm = TRUE), 
    n_years = n_distinct(year),
    .groups = "drop"
  )


######################################################################################################
################################### LOAD DEFAULT QUINCY GPP DATA #####################################
######################################################################################################


setwd("/path/to/quincy/default/")

site_folders_dev <- dir(pattern = "^[A-Z]{2}-[A-Za-z0-9]+")

combined_datasimdev <- list()

for (site in site_folders_dev) {
  
  file_path <- file.path(site, "output", "vegfluxC_daily.txt")
  
  if (file.exists(file_path)) {
    site_data <- read.delim(file_path, sep = "")
    
    site_data$Site <- site
    
    combined_datasimdev[[site]] <- site_data
  } else {
    warning(paste("File not found for site:", site))
  }
}

# Combine all data frames into one big data frame
final_datadev <- do.call(rbind, combined_datasimdev)

# helper vector for assigning years
yearvec_templatedev <- rep(1:124, each = 365)
yearvec_all_sitesdev <- rep(yearvec_templatedev, times = 57)

final_data_cut_dev <- final_datadev %>%
  mutate(year = head(yearvec_all_sitesdev, nrow(.))) %>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(across("Time", round, 0))%>%
  mutate( Date = paste(Time, year + 1900, sep = "/"),  
          Date = as.POSIXct(Date, format = "%j/%Y")       
  ) %>%
  mutate(Site = gsub("-", "_", Site)) 

# merge default GPP data with EC data
matched_data_dev <- final_data_cut_dev %>%
  inner_join(combined_data, by = c("Site" = "site", "Date" = "date"))

complete_years_data_dev_gpp <- matched_data_dev %>%
  group_by(Site, year) %>%
  summarize(
    days_in_year = n(),
    .groups = "drop"
  ) %>%
  filter(days_in_year == 365 | days_in_year == 366) %>%  
  inner_join(matched_data_dev, by = c("Site", "year")) %>%
  select(Site, Time, Date, year, GPP.x, GPP.y) %>%
  mutate(GPP.y = ifelse(GPP.y < 0, 0, GPP.y))

annualsumGPP_dev <- complete_years_data_dev_gpp %>%
  group_by(year, Site) %>%
  summarize(
    total_GPP.x = sum(GPP.x, na.rm = TRUE),
    total_GPP.y = sum(GPP.y, na.rm = TRUE),
    .groups = "drop"  
  ) %>%
  group_by(Site) %>%
  summarize(
    mean_annual_GPP.x = mean(total_GPP.x, na.rm = TRUE),
    sd_annual_GPP.x = sd(total_GPP.x, na.rm = TRUE),  
    mean_annual_GPP.y = mean(total_GPP.y, na.rm = TRUE),
    sd_annual_GPP.y = sd(total_GPP.y, na.rm = TRUE),  
    n_years = n_distinct(year),
    .groups = "drop"
  )


######################################################################################################
################################### LOAD KOPPEN CLIMATE DATA #########################################
######################################################################################################

# raster data from Beck et al. 2018 (doi: 10.1038/sdata.2018.214)
koppen_raster <- rast("/path/to/raster/data/Map_KG-Global/Beck2018_KG/1991_2020/koppen_geiger_0p00833333.tif")

# Create a data frame of coordinates
coords <- data.frame(
  lon = plumber2sites$lon,  
  lat = plumber2sites$lat   
)

# Extract KG classifications 
coords$koppen_class <- extract(koppen_raster, coords[, c("lon", "lat")])[, 2]

# legend converting koppen geiger class numbers to actual names
num2clas<-read.delim("/path/to/legend/koppen_geiger_tif/legend.txt", sep = "", header =F) %>%
  slice(-c(1,2,33:39)) %>%
  select(-c(3:13))%>%
  mutate(across(everything(), ~ gsub(":", "", .)))

colnames(num2clas)[1] <- "koppen_class"
colnames(num2clas)[2] <- "koppen_class_name"
num2clas <- num2clas %>%
  mutate(koppen_class = as.numeric(koppen_class))


# add Köppen-Geiger classifications to plumber 2 site info
plumber2sites_edit <- plumber2sites %>%
  left_join(coords, by = c("lon", "lat")) %>%
  left_join(num2clas, by = "koppen_class"  )

# add site info
plumber2siteschanged <- plumber2sites_edit %>%
  mutate(Site.ID = gsub("-", "_", Site.ID)) %>%
  select(c(1,29,30))


# Add Koeppen-Geiger Classification to GPP data 
annualgpp_siteinfo_sim <- annualsumGPP_sim %>%
  inner_join(plumber2siteschanged, by = c("Site" = "Site.ID"))

# Filter plumber2sites to match sites in annualsumGPP_dev
annualgpp_siteinfo_dev <- annualsumGPP_dev %>%
  inner_join(plumber2siteschanged, by = c("Site" = "Site.ID"))



######################################################################################################
################################### CALCULATE EOS FROM GPP DATA ######################################
######################################################################################################

# Using the package phenofit to extract phenological dates
# from multiple eddy covariance towers

#------------------------------------------
# 1. Set up -----
#------------------------------------------
install.packages("devtools")
library(devtools)
devtools::install_github("geco-bern/FluxDataKit")

# Load packages
packages <- c("tidyverse", "readxl", "here", "zoo", "ggpubr",
              "ncdf4", "lubridate", "FluxDataKit", "phenofit", 
              "gridExtra", "ggrepel")

lapply(packages, library, character.only = TRUE)

# Define the folder containing NetCDF files
nc_folder <- "path/to/PLUMBER/ECflux/files"

#------------------------------------------
# 2. Functions ----
#------------------------------------------

# Prepare flux data (same as before)
prep.flux.data <- function(site, var, start_date) {
  gpp <- ncvar_get(site, var)
  time <- ncvar_get(site, "time")
  start_date <- ymd_hms(start_date)
  dates <- start_date + time
  
  fluxnet_gpp <- data.frame(date = dates, GPP = (gpp * 1e-6 * (60 * 30))) %>%
    mutate(date = as.Date(dates)) %>%
    group_by(date) %>%
    summarise(GPP = sum(GPP, na.rm = TRUE)) %>%
    mutate(year = format(date, "%Y")) %>%
    filter(!(month(date) == 2 & day(date) == 29)) # Remove leap days
  
  return(fluxnet_gpp)
}

get.dates.flux <- function(input_dataset) {
  input <- check_input(input_dataset$date, input_dataset$GPP, south = FALSE, wmin = 0)
  cat("Input data check completed\n")
  print(head(input$y))  # Check the GPP time series input
  
  if (all(is.na(input$y))) {
    warning("All input GPP values are NA!")
    return(NULL)
  }
  
  plot_input(input)
  
  brks_mov <- tryCatch({
    season_mov(input, options = list(
      rFUN = "smooth_wHANTS",
      ntpperyear = "365",
      wFUN = "wTSM",
      iters = 1,
      nf = 2,
      verbose = TRUE
    ))
  }, error = function(e) {
    message("Error in season_mov:", e)
    return(NULL)
  })
  
  if (is.null(brks_mov)) {
    warning("season_mov failed!")
    return(NULL)
  }
  
  plot_season(input, brks_mov, ylab = "GPP (mol C/day)")
  
  fit <- tryCatch({
    curvefits(input, brks_mov, options = list(
      methods = c("Beck", "Elmore"),
      iters = 2,
      wFUN = "wTSM"
    ))
  }, error = function(e) {
    message("Error in curvefits:", e)
    return(NULL)
  })
  
  if (is.null(fit)) {
    warning("curvefits failed!")
    return(NULL)
  }
  
  dfit <- get_fitting(fit)
  g <- plot_curvefits(dfit, brks_mov, NULL, ylab = "GPP", "Time")
  grid::grid.newpage(); grid::grid.draw(g)
  
  pheno <- tryCatch({
    get_pheno(fit, method = c("Beck", "Elmore"), TRS = c(0.1, 0.2, 0.5), IsPlot = TRUE, show.title = FALSE)
  }, error = function(e) {
    message("Error in get_pheno:", e)
    return(NULL)
  })
  
  return(pheno)
}

#------------------------------------------
# 3. Process all files ----
#------------------------------------------

# List all NetCDF files in the folder
nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)

# Initialize an empty list to store results
all_results <- list()

for (nc_file in nc_files) {
  tryCatch({
    nc_data <- nc_open(nc_file)
    site_name <- tools::file_path_sans_ext(basename(nc_file))
    cat("Processing site:", site_name, "\n")
    
    if (!"GPP" %in% names(nc_data$var)) {
      warning(paste("GPP variable not found in file:", nc_file))
      next
    }
    
    extracted_chars <- substr(site_name, 8, 11)
    
    # Print the result
    cat("Starting year:", extracted_chars, "\n")
    
    
    gpp_data <- prep.flux.data(nc_data, "GPP", paste0(extracted_chars,"-01-01 00:00:00"))
    print(head(gpp_data))
    print(dim(gpp_data))
    
    pheno_data <- get.dates.flux(gpp_data)
    if (is.null(pheno_data)) {
      warning(paste("Failed to process phenological dates for site:", site_name))
      next
    }
    
    all_results[[site_name]] <- pheno_data
    nc_close(nc_data)
  }, error = function(e) {
    message("Error processing file:", nc_file, "\n", e)
    next
  })
}


# Combine all results into a list
all_results_combined <- all_results

# get site names
site_idspheno <- plumber2sites$Site.ID

all_results_site_ids <- substr(names(all_results), 1, 6)

# Filter all_results based on matching site IDs
filtered_results_pheno <- all_results[all_results_site_ids %in% site_idspheno]

# simplify list
simplified_results <- lapply(filtered_results_pheno, function(site_data) {
  if (is.list(site_data$doy)) {
    bind_rows(
      site_data$doy$Beck %>% mutate(source = "Beck"),
      site_data$doy$Elmore %>% mutate(source = "Elmore")
    )
  } else {
    stop("Expected 'doy' to contain a list with 'Beck' and 'Elmore' data frames.")
  }
})

# Combine the list into a single long data frame
long_data_frame <- bind_rows(
  lapply(names(simplified_results), function(site_name) {
    
    simplified_results[[site_name]] %>%
      mutate(site = site_name)
  }),
  .id = "list_index" 
)

filtered_long_data_frame <- long_data_frame %>%
  select(site, source, contains(".eos")) %>%
  mutate(site = substr(site, 1, 6))  

long_format_df <- filtered_long_data_frame %>%
  gather(key = "eos_type", value = "eos_value", TRS1.eos, TRS2.eos, TRS5.eos, DER.eos) %>%
  mutate(eos_type = gsub("\\.eos", "", eos_type)) 

summary_df <- long_format_df %>%
  group_by(site, source, eos_type) %>%
  summarise(
    mean_eos = mean(eos_value, na.rm = TRUE),
    sd_eos = sd(eos_value, na.rm = TRUE),
    .groups = "drop"  
  )


##########

# THESE ARE THE FINISHED EOS DATES FROM FLUX DATA 
summary_df <- summary_df %>%
  left_join(plumber2sites_edit[,c(1,4, 30)], by = c("site" = "Site.ID"))

##########


######### EXTRACT DATES FROM GPP DEFAULT/DYNAMIC MODEL #########

# Upload data from folder, individual sites----

cge_new_vegflux_c <-
  read.csv("/path/to/dynamic/simulations/CGE_000/output/vegfluxC_daily.txt",
           sep="")

cge_new_vegdiag <-
  read.csv("/path/to/dynamic/simulations/CGE_000/output/veg_diagnostics_daily.txt",
           sep="")

esmag_new_vegflux_c <-
  read.csv("/path/to/dynamic/simulations/ES-MaG/output/vegfluxC_daily.txt",
           sep="")

esmag_new_vegdiag <-
  read.csv("/path/to/dynamic/simulations/ES-MaG/output/veg_diagnostics_daily.txt",
           sep="")

ustol_new_vegflux_c <-
  read.csv("/path/to/ustol/dynamic/simulations/US-TOL/output/vegfluxC_daily.txt",
           sep="")

ustol_new_vegdiag <-
  read.csv("/path/to/ustol/dynamic/simulations/US-TOL/output/veg_diagnostics_daily.txt",
           sep="")

cge_dev_vegflux_c <-
  read.csv("/path/to/default/simulations/CGE_000/output/vegfluxC_daily.txt",
           sep="")

cge_dev_vegdiag <-
  read.csv("/path/to/default/simulations/CGE_000/output/veg_diagnostics_daily.txt",
           sep="")

esmag_dev_vegflux_c <-
  read.csv("/path/to/default/simulations/ES-MaG/output/vegfluxC_daily.txt",
           sep="")

esmag_dev_vegdiag <-
  read.csv("/path/to/default/simulations/ES-MaG/output/veg_diagnostics_daily.txt",
           sep="")

ustol_dev_vegflux_c <-
  read.csv("/path/to/ustol/default/simulations/US-TOL/output/vegfluxC_daily.txt",
           sep="")

ustol_dev_vegdiag <-
  read.csv("/path/to/ustol/default/simulations/US-TOL/output/veg_diagnostics_daily.txt",
           sep="")

iedri_new_vegflux_c <-
  read.csv("/path/to/dynamic/simulations/IE-Dri/output/vegfluxC_daily.txt",
           sep="")

iedri_new_vegdiag <-
  read.csv("/path/to/dynamic/simulations/IE-Dri/output/veg_diagnostics_daily.txt",
           sep="")

iedri_dev_vegflux_c <-
  read.csv("/path/to/default/simulations/IE-Dri/output/vegfluxC_daily.txt",
           sep="")

iedri_dev_vegdiag <-
  read.csv("/path/to/default/simulations/IE-Dri/output/veg_diagnostics_daily.txt",
           sep="")

sedeg_new_vegflux_c <-
  read.csv("/path/to/dynamic/simulations/SE-Deg/output/vegfluxC_daily.txt",
           sep="")


# these data are from file 'Figure 3 and Table 4.R' and saved as RDS files. They can just be re-used from the environment instead, e.g. usics_flux <- US_ICs_EC_daily
cge_gcc <- 
  readRDS("/path/to/gcc/data/climgrass_gcc.rds")

esmag_flux_gpp <- 
  readRDS("/path/to/gpp/data/ESMaG_ecflux_gpp.rds")

# ustol_flux <- 
#   readRDS("/path/to/gpp/data/ustol_EC_data.rds")

usics_flux <- 
  readRDS("/path/to/gpp/data/US_ICs_EC_daily.rds")

# nc files ----
iedri_flux <-
  nc_open("/path/to/file/IE-Dri_2003-2005_LaThuile_Flux.nc") # object 'plumfluxIEDri' in 'Figure 3 and Table 4.R'

# batch upload ----
base_dir_mod <- 
  "/path/to/dynamic/simulations/"

# List all .dat files in the directory
site_dirs <- 
  list.dirs(base_dir_mod, full.names = TRUE, recursive = FALSE)

all_quincy_results_diag <-
  list()

# Veg diagnostics
# Loop over each site directory
for (site_dir in site_dirs) {
  file_path <- 
    file.path(site_dir, "output", 
              c("veg_diagnostics_daily.txt"))
  
  # Check if the file exists
  if (!file.exists(file_path)) {
    cat("Skipping site:", site_dir, "- file not found.\n")
    next
  }
  
  # Extract the site name from the directory path
  site_name <- basename(site_dir)
  
  cat("Processing site:", site_name, "\n")
  
  # Use tryCatch to handle errors during processing
  tryCatch({
    # Read and process the vegfluxC_daily.txt file
    site_data <- 
      read.delim(file_path, sep = "")
    
    # Save the results in the list with the site name as the key
    all_quincy_results_diag[[site_name]] <- 
      site_data
    
    cat("Successfully processed site:", site_name, "\n")
    
  }, error = function(e) {
    # Handle errors and log them
    cat("Error processing site:", site_name, "\n")
    cat("Error message:", e$message, "\n")
    all_quincy_results_diag[[site_name]] <- list(error = e$message)
  })
}

#------------------------------------------
# Functions ----
#------------------------------------------
# Add dates and convert units in gpp_umol_per_sec to quincy output
prep.quincy.data <- 
  function(dataset) {
    start_date <- 
      ymd("1901-01-01")
    
    num_rows <- 
      nrow(dataset)
    
    date_sequence <- c() # Initialize variables
    current_date <- start_date
    
    while (length(date_sequence) < num_rows) { # Generate date sequence while skipping leap days
      if (!(month(current_date) == 2 & day(current_date) == 29)) {
        date_sequence <- c(date_sequence, current_date)
      }
      current_date <- current_date + days(1)
    }
    
    dataset_final <- 
      dataset %>%
      mutate(date = date_sequence) %>% 
      mutate(date = as.Date(date)) %>% 
      mutate(year = format(date, "%Y"))
    
    return(dataset_final)
  }

prep.flux.data <- 
  function(site, var, start_date) {
    gpp <-
      ncvar_get(site, var)
    
    time <- 
      ncvar_get(site, "time")
    
    start_date <- 
      ymd_hms(start_date)
    
    dates <- 
      start_date + time
    
    fluxnet_gpp <- 
      data.frame(date = dates, GPP = (gpp*1e-6*(60*30))) %>% # Convert in mol/ day
      mutate(date = as.Date(dates)) %>%  # Extract date part
      group_by(date) %>%  # Group by date
      summarise(GPP = sum(GPP, na.rm = TRUE)) %>%
      mutate(year = format(date, "%Y"))
    
    fluxnet_gpp <- 
      fluxnet_gpp %>% # Filter out leap days (February 29)
      filter(!(month(date) == 2 & day(date) == 29))
    
    return(fluxnet_gpp)
  }

get.dates.quincy <- 
  function(input_dataset, option_south) {
    
    input <-
      check_input(input_dataset$date, 
                  input_dataset$GPP, # or GPP or other index (Leaf_C)
                  south = option_south) 
    
    plot_input(input)
    
    brks_mov <- 
      season_mov(input,
                 options = list(
                   rFUN = "smooth_wHANTS", # suggested in Kong 2022 as the most balanced 
                   ntpperyear = "365",
                   wFUN = "wTSM", # suggested in Kong 2022 as the most balanced 
                   iters = 1, # suggested in Kong 2022
                   nf = 2,
                   # r_max = 0.1,
                   # lambda = NULL,
                   # MaxPeaksPerYear = 1,
                   # MaxTroughsPerYear = 2,
                   verbose = TRUE)
                 # years.run = (2003:2015)
      )
    
    plot_season(input, 
                brks_mov, 
                ylab = ("GPP")
    )
    
    fit <- 
      curvefits(input, 
                brks_mov,
                options = list(
                  methods = c("Beck", "Elmore"), #,"Zhang"
                  iters = 2, # default
                  wFUN = "wTSM")
      )
    
    dfit <- get_fitting(fit)
    ## visualization 
    g <- 
      plot_curvefits(dfit, 
                     brks_mov, 
                     NULL, 
                     ylab = "GPP", 
                     "Time")
    
    grid::grid.newpage(); grid::grid.draw(g) # plot to check the curve fitting
    
    pheno <- 
      get_pheno(fit, 
                method = c("Beck","Elmore"), # recommended in Kong 2022
                TRS = c(0.1,0.2,0.5, 0.6), 
                IsPlot = TRUE, 
                show.title = FALSE)
    
    return(pheno)
    
  }

get.dates.gcc <- 
  function(input_dataset, option_south) {
    
    input <-
      check_input(input_dataset$date, 
                  input_dataset$gcc_90, # or GPP or other index (Leaf_C)
                  south = option_south) 
    
    plot_input(input)
    
    brks_mov <- 
      season_mov(input,
                 options = list(
                   rFUN = "smooth_wHANTS", # suggested in Kong 2022 as the most balanced 
                   ntpperyear = "365",
                   wFUN = "wTSM", # suggested in Kong 2022 as the most balanced 
                   iters = 1, # suggested in Kong 2022
                   nf = 2,
                   # r_max = 0.1,
                   # lambda = NULL,
                   # MaxPeaksPerYear = 1,
                   # MaxTroughsPerYear = 2,
                   verbose = TRUE)
                 # years.run = (2003:2015)
      )
    
    plot_season(input, 
                brks_mov, 
                ylab = ("GPP")
    )
    
    fit <- 
      curvefits(input, 
                brks_mov,
                options = list(
                  methods = c("Beck", "Elmore"), #,"Zhang"
                  iters = 2, # default
                  wFUN = "wTSM")
      )
    
    dfit <- get_fitting(fit)
    ## visualization 
    g <- 
      plot_curvefits(dfit, 
                     brks_mov, 
                     NULL, 
                     ylab = "GPP", 
                     "Time")
    
    grid::grid.newpage(); grid::grid.draw(g) # plot to check the curve fitting
    
    pheno <- 
      get_pheno(fit, 
                method = c("Beck","Elmore"), # recommended in Kong 2022
                TRS = c(0.1,0.2,0.5, 0.6), 
                IsPlot = TRUE, 
                show.title = FALSE)
    
    return(pheno)
    
  }

get.dates.quincy.out <- 
  function(input_dataset) {
    
    # Identify years to exclude (where growing_season is all 1s or all 0s)
    excluded_years <- input_dataset %>%
      group_by(year) %>%
      filter(all(growing_season == 1) | all(growing_season == 0)) %>%
      summarise(excluded_reason = ifelse(all(growing_season == 1), "All 1s", "All 0s")) %>%
      ungroup()
    
    # Assign a unique season ID to each continuous period where growing_season == 1
    input_dataset <- input_dataset %>%
      arrange(date) %>%  # Ensure data is ordered chronologically
      mutate(season_id = cumsum(growing_season == 1 & lag(growing_season, default = 0) == 0))
    
    # Find SOS and EOS for each growing season
    pheno_quincy_output <- input_dataset %>%
      group_by(season_id) %>%
      summarise(
        sos_quincy_out = min(date[growing_season == 1]),
        eos_quincy_out = max(date[growing_season == 1])
      ) %>%
      ungroup() %>% 
      mutate(year = as.character(lubridate::year(sos_quincy_out))) %>%  # Extract year from SOS
      select(year, sos_quincy_out, eos_quincy_out)   
    
    return(list(
      pheno_quincy_output = pheno_quincy_output,
      excluded_years = excluded_years
    ))
  }

compare.dates.flux <- 
  function(pheno, pheno_quincy, pheno_quincy_output) { 
    
    pheno_beck_flux <-
      pheno[["date"]][["Beck"]] %>% 
      rename_with(~ paste0(., "_flux_beck"), -c(flag, origin))
    
    pheno_elmore_flux <-
      pheno[["date"]][["Elmore"]] %>% 
      rename_with(~ paste0(., "_flux_elmore"), -c(flag, origin))
    
    pheno_elmore_quincy <-
      pheno_quincy[["date"]][["Elmore"]] %>% 
      rename_with(~ paste0(., "_quincy_elmore"), -c(flag, origin))
    
    pheno_beck_quincy <-
      pheno_quincy[["date"]][["Beck"]] %>% 
      rename_with(~ paste0(., "_quincy_beck"), -c(flag, origin))
    
    pheno_compare_beck <-
      full_join(pheno_beck_flux,
                pheno_beck_quincy)
    
    pheno_compare_elmore <-
      full_join(pheno_elmore_flux,
                pheno_elmore_quincy) 
    
    pheno_compare_phenofit <-
      full_join(pheno_compare_beck,
                pheno_compare_elmore) %>% 
      mutate(year = format(origin, "%Y"))
    
    pheno_compare_full <-
      full_join(pheno_compare_phenofit,
                pheno_quincy_output)
    
    return(pheno_compare_full)
  }


#------------------------------------------
# Edit data ----
#------------------------------------------
# individual sites----
# Add the datetime column to existing data frame
# convert data into umol/sec
sites <- 
  list(
    cge_quincy_c = cge_new_vegflux_c,
    cge_quincy_diag = cge_new_vegdiag,
    esmag_quincy_c = esmag_new_vegflux_c,
    esmag_quincy_diag = esmag_new_vegdiag,
    ustol_quincy_c = ustol_new_vegflux_c,
    ustol_quincy_diag = ustol_new_vegdiag,
    iedri_quincy_c = iedri_new_vegflux_c,
    cge_dev_quincy_c = cge_dev_vegflux_c,
    cge_dev_quincy_diag = cge_dev_vegdiag,
    esmag_dev_quincy_c = esmag_dev_vegflux_c,
    esmag_dev_quincy_diag = esmag_dev_vegdiag,
    ustol_dev_quincy_c = ustol_dev_vegflux_c,
    ustol_dev_quincy_diag = ustol_dev_vegdiag,
    iedri_dev_quincy_c = iedri_dev_vegflux_c
  )

# Apply prep.quincy.data to each element in the list
quincy_sites <- 
  lapply(sites, 
         prep.quincy.data)

# Assign the processed data back to the original variable names
list2env(quincy_sites, 
         envir = .GlobalEnv)

# add year col tuo flux dataset
esmag_flux_gpp <-
  esmag_flux_gpp %>% 
  mutate(year = format(date, "%Y"))

usics_flux_gpp <-
  usics_flux %>% 
  mutate(year = format(Date, "%Y"),
         GPP = abs(GPP_NT),
         date = as.Date(Date)) %>% 
  filter(!is.na(Date))

# process nc data
iedri_flux_gpp <-
  prep.flux.data (iedri_flux,
                  "GPP",
                  "2003-01-01 00:00:00")

# get sos/eos
# cge
# methods
pheno_quincy_cge <-
  get.dates.quincy(cge_quincy_c, FALSE)

# quincy output
pheno_quincy_out_cge <-
  get.dates.quincy.out(cge_quincy_diag)

# gcc
pheno_gcc_cge <-
  get.dates.gcc(cge_gcc, FALSE)

# compare
pheno_compare_cge <-
  compare.dates.flux( 
    pheno_gcc_cge,
    pheno_quincy_cge,
    pheno_quincy_out_cge$pheno_quincy_output)

# us tol
# methods
pheno_quincy_ustol <-
  get.dates.quincy(ustol_quincy_c, FALSE)

# # quincy output
pheno_quincy_out_ustol <-
  get.dates.quincy.out(ustol_quincy_diag)

# flux
pheno_flux_usics <-
  get.dates.quincy(usics_flux_gpp, FALSE)

# compare
pheno_compare_ustol <-
  compare.dates.flux(pheno_flux_usics,
                     pheno_quincy_ustol, 
                     pheno_quincy_out_ustol$pheno_quincy_output)

# es mag
# methods
pheno_quincy_esmag <-
  all_quincy_results[["ES-MaG"]]

pheno_quincy_esmag <-
  get.dates.quincy(esmag_quincy_c, FALSE)

# quincy output
pheno_quincy_out_esmag <-
  dates_grass_out_new[["ES-MaG"]]

# flux
pheno_flux_esmag <-
  get.dates.quincy(esmag_flux_gpp, FALSE)

# compare
pheno_compare_esmag <-
  compare.dates.flux( 
    pheno_flux_esmag,
    pheno_quincy_esmag,
    pheno_quincy_out_esmag$pheno_quincy_output)

# cge dev
# methods
pheno_quincy_cge_dev <-
  get.dates.quincy(cge_dev_quincy_c, FALSE)

# quincy output
pheno_quincy_out_cge_dev <-
  get.dates.quincy.out(cge_dev_quincy_diag)

# compare
pheno_compare_cge_dev <-
  compare.dates.flux( 
    pheno_gcc_cge,
    pheno_quincy_cge_dev,
    pheno_quincy_out_cge_dev$pheno_quincy_output)

# us tol dev
# methods
pheno_quincy_ustol_dev <-
  get.dates.quincy(ustol_dev_quincy_c, FALSE)

# quincy output
pheno_quincy_out_ustol_dev <-
  get.dates.quincy.out(ustol_dev_quincy_diag)

# compare
pheno_compare_ustol_dev <-
  compare.dates.flux(pheno_flux_usics,
                     pheno_quincy_ustol_dev, 
                     pheno_quincy_out_ustol_dev$pheno_quincy_output)

# es mag dev
# methods
pheno_quincy_esmag_dev <-
  get.dates.quincy(esmag_dev_quincy_c, FALSE)

# quincy output
pheno_quincy_out_esmag_dev <-
  get.dates.quincy.out(esmag_dev_quincy_diag)

# flux
pheno_flux_esmag <-
  get.dates.quincy(esmag_flux_gpp, FALSE)

# compare
pheno_compare_esmag_dev <-
  compare.dates.flux( 
    pheno_flux_esmag,
    pheno_quincy_esmag_dev,
    pheno_quincy_out_esmag_dev$pheno_quincy_output)

# ie dri
# flux
pheno_flux_iedri <-
  get.dates.quincy(iedri_flux_gpp, FALSE)

# compare
pheno_compare_iedri <-
  compare.dates.flux(pheno_flux_iedri,
                     pheno_grass_new$`IE-Dri`, 
                     dates_grass_out_new$`IE-Dri`$pheno_quincy_output)

# ie dri dev
# quincy output
pheno_quincy_out_iedri_dev <-
  get.dates.quincy.out(iedri_dev_quincy_diag)

# compare
pheno_compare_iedri_dev <-
  compare.dates.flux(pheno_flux_iedri,
                     dates_grass_dev$`IE-Dri`, 
                     pheno_quincy_out_iedri_dev$pheno_quincy_output)

# batch edit----
# dates from quincy output
all_quincy_results_diag <-
  all_quincy_results_diag %>% 
  lapply(., prep.quincy.data)

dates_quincy_out <-
  all_quincy_results_diag %>% 
  lapply(., get.dates.quincy.out)



# get dates from from GPP

# redo this for default model as well as dynamic model. 
# dynamic model: all_quincy_results_mod
# default model: all_quincy_results


# to get whether sites are in the northern or southern hemisphere
# List all .nc files in the directory
all_quincy_results <-
  list()

nc_files <- 
  list.files(
    path = "path/to/PLUMBER/flux/data", 
    pattern = "\\.nc$", 
    full.names = TRUE
  )

# Extract site, start date, and end date
nc_data <- 
  data.frame(
    file_path = nc_files,
    stringsAsFactors = FALSE
  ) %>%
  mutate(
    site = sub("_.*", "", basename(file_path)),                
    date_range = sub(".*_(\\d{4}-\\d{4})_.*", "\\1", file_path),
    start_date = as.numeric(sub("-.*", "", date_range)),    
    end_date = as.numeric(sub(".*-", "", date_range))       
  ) %>%
  select(site, start_date, end_date) 


# List all site dirs
site_dirs <- 
  list.dirs(base_dir_mod, full.names = TRUE, recursive = FALSE)

ec_sites_temp <- 
  basename(site_dirs)  

# Filter nc_data to include only matching sites
nc_data <- 
  nc_data %>%
  filter(site %in% ec_sites_temp)


plumber2coord <- 
  read.delim(file = "/read/site/info/latlon/all_sites_list_t.dat", 
             sep = "") %>%
  select(c(1,3)) %>%
  mutate(site = Site.ID,
         south = ifelse(lat < 0, TRUE, FALSE)) %>%
  select(-c(Site.ID, lat))

nc_data <- 
  left_join(nc_data, 
            plumber2coord %>% 
              select(site, south), by = "site")

# compute dates
for (site_dir in site_dirs) {
  # Define path
  vegflux_path <- 
    file.path(site_dir, "output", "vegfluxC_daily.txt")
  
  if (!file.exists(vegflux_path)) {
    cat("Skipping site:", site_dir, "- vegfluxC_daily.txt not found.\n")
    next
  }
  
  # get site name
  site_name <- basename(site_dir)
  
  # Filter nc_data to get the start and end years for current site
  site_nc_data <- nc_data %>% filter(site == site_name)
  if (nrow(site_nc_data) == 0) {
    cat("No matching NC data found for site:", site_name, "\n")
    next
  }
  
  start_year <- as.numeric(site_nc_data$start_date)
  end_year <- as.numeric(site_nc_data$end_date)
  
  cat("Processing site:", site_name, "\n")
  
  tryCatch({
    site_data <- 
      read.delim(vegflux_path, sep = "") 
    prepped_data <- 
      prep.quincy.data(site_data)
    
    prepped_data <- prepped_data %>% 
      filter(year >= start_year & year <= end_year)
    
    pheno_data <- 
      get.dates.quincy(prepped_data , 
                       option_south = nc_data$south)
    
    # Save results
    all_quincy_results[[site_name]] <- 
      # pheno_data
      
      cat("processed site:", site_name, "\n")
    
  }, error = function(e) {
    # Handle errors and log them
    cat("Error processing site:", site_name, "\n")
    cat("Error message:", e$message, "\n")
    all_quincy_results[[site_name]] <- list(error = e$message)
  })
}


# default model EOS results: 
all_quincy_results

# dynamic model EOS results:
all_quincy_results_mod


#-----------------------------

# Combine the default/dynamic data 'all_quincy_results'/'all_quincy_results_mod' with the EC flux data 'summary_df'

#--------------------

# dynamic model data:

#-----------------------------

# get site names
quincy_mod_site_ids <- substr(names(all_quincy_results_mod), 1, 6)

# Filter all_results based on matching site IDs
filtered_results_quincy_mod <- all_quincy_results_mod[quincy_mod_site_ids %in% site_idspheno]

simplified_results_quincy_mod <- lapply(filtered_results_quincy_mod, function(site_data) {
  if (is.list(site_data$doy)) {
    return(site_data$doy)
  } else {
    stop("Expected 'doy'")
  }
})

# Change to long df from list
long_data_frame_qm <- bind_rows(
  lapply(names(simplified_results_quincy_mod), function(site_name) {
    site_data <- simplified_results_quincy_mod[[site_name]]
  }),
  .id = "list_index"
)


# Filter for EOS only
filtered_long_data_frame_qm <- long_data_frame_qm %>%
  select(site, source, contains(".eos")) %>%
  mutate(site = substr(site, 1, 6))  

# select relevant columns
long_format_df_qm <- filtered_long_data_frame_qm %>%
  select(c("site","source","DER.eos")) 

# calc. mean and sd 
summary_df_qm <- long_format_df_qm %>%
  group_by(site) %>%
  summarise(
    mean_eos = mean(DER.eos, na.rm = TRUE),
    sd_eos = sd(DER.eos, na.rm = TRUE),
    .groups = "drop"
  )

# add site info
summary_df_qm <- summary_df_qm %>%
  left_join(plumber2sites_edit[,c(1,4, 30)], by = c("site" = "Site.ID"))

# filter EC data for EOS
summary_df_DER <- summary_df %>%
  filter(metric == 'DER.eos')

# Join with flux data
merged_df_fin <- summary_df_qm %>%
  inner_join(summary_df_DER, by = c("site", "PFT", "koppen_class_name"), 
             suffix = c(".qm", ".ec"))



#-----------------------------

# repeat edit for default model


#-----------------------------
# get site name
quincy_dev_site_ids <- substr(names(all_quincy_results), 1, 6)

# Filter all_results based on matching site IDs
filtered_results_quincy_dev <- all_quincy_results[quincy_dev_site_ids %in% site_idspheno]

simplified_results_quincy_dev <- lapply(filtered_results_quincy_dev, function(site_data) {
  if (is.list(site_data$doy)) {
    return(site_data$doy)
  } else {
    stop("Expected 'doy' to contain a list with 'Beck' and 'Elmore' components.")
  }
})

# Check the structure of the simplified list
str(simplified_results_quincy_dev)

# Reformat to long df from list
long_data_frame_dev <- bind_rows(
  lapply(names(simplified_results_quincy_dev), function(site_name) {
    site_data <- simplified_results_quincy_dev[[site_name]]
  }),
  .id = "list_index"
)

# Filter for EOS
filtered_long_data_frame_dev <- long_data_frame_dev %>%
  select(site, source, contains(".eos")) %>%
  mutate(site = substr(site, 1, 6))

# select relevant columns
long_format_df_dev <- filtered_long_data_frame_dev %>%  
  select(c("site","source","DER.eos")) 

# calc. mean and sd
summary_df_dev <- long_format_df_dev %>%
  group_by(site) %>%
  summarise(
    mean_eos = mean(DER.eos, na.rm = TRUE),
    sd_eos = sd(DER.eos, na.rm = TRUE),
    .groups = "drop" 
  )

# add site info
summary_df_dev <- summary_df_dev %>%
  left_join(plumber2sites_edit[,c(1,4, 30)], by = c("site" = "Site.ID"))

# join with EC EOS data
merged_df_fin_dev <- summary_df_dev %>%
  inner_join(summary_df_DER, by = c("site", "PFT", "koppen_class_name"), 
             suffix = c(".dev", ".ec"))



######################################################################################################
####################################### CREATE FIGURE 5 ##############################################
######################################################################################################

# Figure 5a
developmulti <- ggplot() +
  
  geom_errorbar(
    data=annualgpp_siteinfo_dev%>% filter(!is.na(sd_annual_GPP.x) & !is.na(sd_annual_GPP.y)),
    aes(
      x = mean_annual_GPP.x*12.0107, 
      ymin = mean_annual_GPP.y*12.0107 - sd_annual_GPP.y*12.0107, 
      ymax = mean_annual_GPP.y*12.0107 + sd_annual_GPP.y*12.0107, 
      color = koppen_class_name 
    ), 
    width = 45, alpha = 0.8, lwd = 1
  ) +
  geom_errorbarh(
    data=annualgpp_siteinfo_dev %>% filter(!is.na(sd_annual_GPP.x) & !is.na(sd_annual_GPP.y)),
    aes(
      y = mean_annual_GPP.y*12.0107,
      xmin = mean_annual_GPP.x*12.0107 - sd_annual_GPP.x*12.0107, 
      xmax = mean_annual_GPP.x*12.0107 + sd_annual_GPP.x*12.0107,
      color = koppen_class_name 
    ), height = 45, alpha = 0.8, lwd = 1) +
  
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  
  geom_point(data=annualgpp_siteinfo_dev,
             aes(
               x = (mean_annual_GPP.x*12.0107), 
               y = (mean_annual_GPP.y*12.0107),
               color = koppen_class_name, 
               shape = koppen_class_name), size = 8, stroke = 2) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_text(size = 36, face = 'bold'),
    axis.title.y = element_markdown(size = 36, face = 'plain'),
    axis.title.x = element_markdown(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   
    axis.text.y = element_text(size = 36, face = 'plain'),   
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) +
  
  scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, 500))+
  scale_x_continuous(limits = c(0, 3000), breaks = seq(0, 3000, 500))+
  coord_cartesian(xlim = c(0, 3000), ylim = c(0, 3000))+
  
  scale_color_discrete(name = "Köppen-Geiger Classification", type = pals::glasbey())+ 
  scale_shape_manual(name = "Köppen-Geiger Classification",values = shapes) +
  
  guides(
    color = guide_legend(override.aes = list(size = 5), nrow =2, byrow=T),  
    shape = guide_legend(override.aes = list(size = 5), nrow =2, byrow=T)   
  ) +
  
  labs(
    title = "Default Model",
    x = "QUINCY mean annual GPP g C m<sup>-2</sup> y<sup>-1</sup>",
    y = "EC flux mean annual GPP g C m<sup>-2</sup> y<sup>-1</sup>"); print(developmulti)


# Figure 5b

newmod <- ggplot() +
  
  geom_errorbar(data=annualgpp_siteinfo_sim%>% filter(!is.na(sd_annual_GPP.x) & !is.na(sd_annual_GPP.y)),
                aes(
                  x = mean_annual_GPP.x * 12.0107,
                  ymin = mean_annual_GPP.y * 12.0107 - sd_annual_GPP.y * 12.0107,
                  ymax = mean_annual_GPP.y * 12.0107 + sd_annual_GPP.y * 12.0107,
                  color = koppen_class_name
                ),
                width = 45, alpha = 0.8, lwd = 1
  ) +
  
  geom_errorbarh(data=annualgpp_siteinfo_sim%>% filter(!is.na(sd_annual_GPP.x) & !is.na(sd_annual_GPP.y)),
                 aes(
                   y = mean_annual_GPP.y * 12.0107,
                   xmin = mean_annual_GPP.x * 12.0107 - sd_annual_GPP.x * 12.0107,
                   xmax = mean_annual_GPP.x * 12.0107 + sd_annual_GPP.x * 12.0107,
                   color = koppen_class_name
                 ),
                 height = 45, alpha = 0.8, lwd = 1
  ) +
  
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  
  geom_point(data=annualgpp_siteinfo_sim%>% filter(!is.na(sd_annual_GPP.x) & !is.na(sd_annual_GPP.y)),
             aes(
               x = mean_annual_GPP.x * 12.0107,
               y = mean_annual_GPP.y * 12.0107,
               color = koppen_class_name,
               shape = koppen_class_name
             ),
             size = 8, stroke = 2
  ) +
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_text(size = 36, face = 'bold'),
    axis.title.y = element_markdown(size = 36, face = 'plain'),
    axis.title.x = element_markdown(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   
    axis.text.y = element_text(size = 36, face = 'plain'),   
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) +
  
  scale_y_continuous(limits = c(-1000, 3000), breaks = seq(0, 3000, 500)) +
  scale_x_continuous(limits = c(-1000, 3000), breaks = seq(0, 3000, 500)) +
  coord_cartesian(xlim = c(0, 3000), ylim = c(0, 3000))+
  
  scale_color_discrete(
    name = "Köppen-Geiger Classification",
    type = pals::glasbey()
  ) +
  
  scale_shape_manual(
    name = "Köppen-Geiger Classification",
    values = shapes
  ) +
  
  guides(
    color = guide_legend(override.aes = list(size = 5), nrow =2, byrow=T),
    shape = guide_legend(override.aes = list(size = 5), nrow =2, byrow=T)
  ) +
  
  labs(
    title = "Dynamic Model",
    x = "QUINCY mean annual GPP g C m<sup>-2</sup> y<sup>-1</sup>",
    y = "EC flux mean annual GPP g C m<sup>-2</sup> y<sup>-1</sup>"
  );print(newmod)


# Figure 5c
devsen <- ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_errorbar(data = merged_df_fin_dev %>% filter(!site == 'SD-Dem') %>% filter(!is.na(sd_eos.ec) & !is.na(sd_eos.dev)), 
                aes(x = mean_eos.dev, ymin = mean_eos.ec - sd_eos.ec, ymax = mean_eos.ec + sd_eos.ec, 
                    color = as.factor(koppen_class_name)), width = 4, alpha = 0.8, lwd = 1) +
  geom_errorbarh(data = merged_df_fin_dev %>% filter(!site == 'SD-Dem') %>% filter(!is.na(sd_eos.ec) & !is.na(sd_eos.dev)),  
                 aes(y = mean_eos.ec, xmin = mean_eos.dev - sd_eos.dev, xmax = mean_eos.dev + sd_eos.dev, 
                     color = as.factor(koppen_class_name)), height = 4, alpha = 0.8, lwd = 1) +
  geom_point(data = merged_df_fin_dev %>% filter(!site == 'SD-Dem'), aes(x = mean_eos.dev, y = mean_eos.ec, color = as.factor(koppen_class_name), shape = as.factor(koppen_class_name)), 
             size = 8, stroke = 2) +
  labs(
    title = "Default Model",
    x = "QUINCY mean EOS DOY",
    y = "Flux mean EOS DOY"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_text(size = 36, face = 'bold'),
    axis.title.y = element_text(size = 36, face = 'plain'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   
    axis.text.y = element_text(size = 36, face = 'plain'),   
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) +
  scale_y_continuous(limits = c(-50, 500), breaks = c(0, 80, 160, 240, 320, 365)) +
  scale_x_continuous(limits = c(-50, 500), breaks = c(0, 80, 160, 240, 320, 365)) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 365))+
  scale_color_discrete(name = "Köppen-Geiger Classification", type = pals::glasbey())+ 
  scale_shape_manual(name = "Köppen-Geiger Classification",values = shapes) +
  guides(
    color = guide_legend(override.aes = list(size = 5), nrow =2, byrow=T),  
    shape = guide_legend(override.aes = list(size = 5), nrow =2, byrow=T)   
  ); print(devsen)  


# Figure 5d

simsen <- ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_errorbar(data = merged_df_fin %>% filter(!site == 'SD-Dem') %>% filter(!is.na(sd_eos.ec) & !is.na(sd_eos.qm)), 
                aes(x = mean_eos.qm, ymin = mean_eos.ec - sd_eos.ec, ymax = mean_eos.ec + sd_eos.ec, 
                    color = as.factor(koppen_class_name)), width = 4, alpha = 0.8, lwd = 1) +
  geom_errorbarh(data = merged_df_fin %>% filter(!site == 'SD-Dem') %>% filter(!is.na(sd_eos.ec) & !is.na(sd_eos.qm)), 
                 aes(y = mean_eos.ec, xmin = mean_eos.qm - sd_eos.qm, xmax = mean_eos.qm + sd_eos.qm, 
                     color = as.factor(koppen_class_name)), height = 4, alpha = 0.8, lwd = 1) +
  geom_point(data = merged_df_fin%>% filter(!site == 'SD-Dem'), 
             aes(x = mean_eos.qm, y = mean_eos.ec, color = as.factor(koppen_class_name), shape = as.factor(koppen_class_name)), 
             size = 8, stroke = 2) +
  labs(
    title = "Dynamic Model",
    x = "QUINCY mean EOS DOY",
    y = "Flux mean EOS DOY"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA),
    plot.title = element_text(size = 36, face = 'bold'),
    axis.title.y = element_text(size = 36, face = 'plain'),
    axis.title.x = element_text(size = 36, face = 'plain'),
    axis.text.x = element_text(size = 36, face = 'plain'),   
    axis.text.y = element_text(size = 36, face = 'plain'),   
    legend.text = element_text(size = 36, face = 'plain'),
    legend.title = element_text(size = 36, face = 'bold')
  ) +
  scale_y_continuous(limits = c(0, 400), breaks = c(0, 80, 160, 240, 320, 365)) +
  scale_x_continuous(limits = c(0, 400), breaks = c(0, 80, 160, 240, 320, 365)) +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 365))+
  scale_color_discrete(name = "Köppen-Geiger Classification", type = pals::glasbey()) + 
  scale_shape_manual(name = "Köppen-Geiger Classification",values = shapes) +
  guides(
    color = guide_legend(override.aes = list(size = 5, stroke = 1), nrow =2, byrow=T), 
    shape = guide_legend(override.aes = list(size = 5, stroke = 1), nrow =2, byrow=T)  
  ) ; print(simsen)  



# Arrange the plots
combined_plot <- ggarrange(devsen, simsen, 
                           common.legend = TRUE, 
                           legend = "right", 
                           labels = c("a)", "b)", 'c)'),  
                           label.x = 0,
                           label.y = 1,
                           font.label = list(size = 36)); print(combined_plot)



combined_plot2 <- ggarrange(developmulti, NULL, newmod, 
                            NULL, NULL, NULL,
                            devsen, NULL, simsen,
                            nrow = 4, ncol = 3,
                            widths = c(1, .2, 1),
                            heights = c(1,.2, 1, .1),
                            common.legend = TRUE, 
                            legend = "bottom", 
                            labels = c("a)",'', "b)",'','','', 'c)','','d)'),  
                            label.x = 0,
                            label.y = 1,
                            font.label = list(size = 36)); print(combined_plot2)


ggsave(
  "/path/to/figure.png",
  plot = combined_plot2,
  scale = 4,
  width = 2200,
  height = 2000,
  units = "px")






######################################################################################################
####################################### STATS TABLE 5 AND TABLE C2 ###################################
######################################################################################################

r_squared <- function(observed, predicted) {
  model = lm(observed~predicted) 
  summary(model)
  summary(model)$r.squared 
}

library(hydroGOF) 
library(INDperform)


##############

####### Table 5 #######

##############

# ROOT MEAN SQUARE ERROR FOR ALL SITES DEFAULT AND DYNAMIC MODEL


merged_df_fin_fast_sen <- merged_df_fin %>%
  filter(eos_type == "DER")

hydroGOF::rmse(obs = merged_df_fin_fast_sen$mean_eos.ec, sim = merged_df_fin_fast_sen$mean_eos.qm)

hydroGOF::mae(obs = merged_df_fin_fast_sen$mean_eos.ec, sim = merged_df_fin_fast_sen$mean_eos.qm)

INDperform::nrmse(obs = merged_df_fin_fast_sen$mean_eos.ec, pred = merged_df_fin_fast_sen$mean_eos.qm)


merged_df_fin_dev_fast_sen <-  merged_df_fin_dev %>%
  filter(eos_type == "DER")

hydroGOF::rmse(obs = merged_df_fin_dev_fast_sen$mean_eos.ec, sim = merged_df_fin_dev_fast_sen$mean_eos.dev)

hydroGOF::mae(obs = merged_df_fin_dev_fast_sen$mean_eos.ec, sim = merged_df_fin_dev_fast_sen$mean_eos.dev)

INDperform::nrmse(obs = merged_df_fin_dev_fast_sen$mean_eos.ec, pred = merged_df_fin_dev_fast_sen$mean_eos.dev)



##############

####### Table C2 #######

##############



# MAE, RMSE and NRMSE

# new turnover

hydroGOF::rmse(obs = annualgpp_siteinfo_sim$mean_annual_GPP.y, sim = annualgpp_siteinfo_sim$mean_annual_GPP.x)

hydroGOF::mae(obs = annualgpp_siteinfo_sim$mean_annual_GPP.y, sim = annualgpp_siteinfo_sim$mean_annual_GPP.x)

INDperform::nrmse(obs = annualgpp_siteinfo_sim$mean_annual_GPP.y, pred = annualgpp_siteinfo_sim$mean_annual_GPP.x)


# default turnover

hydroGOF::mae(obs = annualgpp_siteinfo_dev$mean_annual_GPP.y , sim = annualgpp_siteinfo_dev$mean_annual_GPP.x)

hydroGOF::rmse(obs = annualgpp_siteinfo_dev$mean_annual_GPP.y , sim = annualgpp_siteinfo_dev$mean_annual_GPP.x)

INDperform::nrmse(obs = annualgpp_siteinfo_dev$mean_annual_GPP.y , pred = annualgpp_siteinfo_dev$mean_annual_GPP.x)


# pearson rho

# new turnover

cor(annualgpp_siteinfo_sim$mean_annual_GPP.y , annualgpp_siteinfo_sim$mean_annual_GPP.x)

# default turnover

cor(annualgpp_siteinfo_dev$mean_annual_GPP.y , annualgpp_siteinfo_dev$mean_annual_GPP.x)



# new turnover

r_squared(observed = annualgpp_siteinfo_sim$mean_annual_GPP.y, predicted = annualgpp_siteinfo_sim$mean_annual_GPP.x)

lmnew <- lm(annualgpp_siteinfo_sim$mean_annual_GPP.y ~ annualgpp_siteinfo_sim$mean_annual_GPP.x)
summary(lmnew)

# default turnover

r_squared(observed = annualgpp_siteinfo_dev$mean_annual_GPP.y, predicted =  annualgpp_siteinfo_dev$mean_annual_GPP.x)

lmdev <- lm(annualgpp_siteinfo_dev$mean_annual_GPP.y ~ annualgpp_siteinfo_dev$mean_annual_GPP.x)
summary(lmdev)


annualgpp_siteinfo_sim_gramm <- annualgpp_siteinfo_sim %>%
  rename(gpp_ec = mean_annual_GPP.y,
         gpp_qm = mean_annual_GPP.x) %>%
  mutate(gpp_ec = gpp_ec * 12.0107,
         gpp_qm = gpp_qm * 12.0107)

annualgpp_siteinfo_dev_gramm <- annualgpp_siteinfo_dev %>%
  rename(gpp_ec = mean_annual_GPP.y,
         gpp_dev = mean_annual_GPP.x) %>%
  mutate(gpp_ec = gpp_ec * 12.0107,
         gpp_dev = gpp_dev * 12.0107)


overall_gpp_stats_for_kop_gramm <- data.frame(
  koppen_class_name = "All Sites",
  n_sites_per_class = NA,
  n_site_years_per_class = sum(annualgpp_siteinfo_sim_gramm$n_years, na.rm = TRUE),
  mean_gpp_ec  = mean(annualgpp_siteinfo_sim_gramm$gpp_ec, na.rm = TRUE),
  mean_gpp_dev = mean(annualgpp_siteinfo_dev_gramm$gpp_dev, na.rm = TRUE),
  mean_gpp_qm  = mean(annualgpp_siteinfo_sim_gramm$gpp_qm, na.rm = TRUE),
  
  MAE.dev   = hydroGOF::mae(obs = annualgpp_siteinfo_dev_gramm$gpp_ec, sim = annualgpp_siteinfo_dev_gramm$gpp_dev),
  RMSE.dev  = hydroGOF::rmse(obs = annualgpp_siteinfo_dev_gramm$gpp_ec, sim = annualgpp_siteinfo_dev_gramm$gpp_dev),
  NRMSE.dev = hydroGOF::nrmse(obs = annualgpp_siteinfo_dev_gramm$gpp_ec, sim = annualgpp_siteinfo_dev_gramm$gpp_dev, norm = 'maxmin'),
  
  MAE.new   = hydroGOF::mae(obs = annualgpp_siteinfo_sim_gramm$gpp_ec, sim = annualgpp_siteinfo_sim_gramm$gpp_qm),
  RMSE.new  = hydroGOF::rmse(obs = annualgpp_siteinfo_sim_gramm$gpp_ec, sim = annualgpp_siteinfo_sim_gramm$gpp_qm),
  NRMSE.new = hydroGOF::nrmse(obs = annualgpp_siteinfo_sim_gramm$gpp_ec, sim = annualgpp_siteinfo_sim_gramm$gpp_qm, norm = 'maxmin')
) 

rmse_new_kopp_class_gpp <- annualgpp_siteinfo_sim_gramm %>%
  group_by(koppen_class_name) %>%
  summarize(mean_gpp_qm = mean(gpp_qm, na.rm=T),
            mean_gpp_ec = mean(gpp_ec, na.rm=T),
            MAE   = hydroGOF::mae(obs = gpp_ec, sim = gpp_qm),
            RMSE  = hydroGOF::rmse(obs = gpp_ec, sim = gpp_qm),
            NRMSE = hydroGOF::nrmse(obs = gpp_ec, sim = gpp_qm, norm = 'maxmin'),
            n_sites_per_class = n(),
            n_site_years_per_class = sum(n_years, na.rm = TRUE)
  ) 


rmse_dev_kopp_class_gpp <- annualgpp_siteinfo_dev_gramm %>%
  group_by(koppen_class_name) %>% 
  summarize(mean_gpp_dev = mean(gpp_dev, na.rm=T),
            mean_gpp_ec = mean(gpp_ec, na.rm=T),
            MAE   = hydroGOF::mae(obs = gpp_ec, sim = gpp_dev),
            RMSE  = hydroGOF::rmse(obs = gpp_ec, sim = gpp_dev),
            NRMSE = hydroGOF::nrmse(obs = gpp_ec, sim = gpp_dev, norm = 'maxmin'),
            n_sites_per_class = n(),
            n_site_years_per_class = sum(n_years, na.rm = TRUE)
  )



# Step 1: Join dev and new model RMSE/MAE tables
rmse_total_by_kop_gpp <- rmse_dev_kopp_class_gpp %>%
  inner_join(
    rmse_new_kopp_class_gpp,
    by = c("koppen_class_name", "n_sites_per_class", "mean_gpp_ec"),
    suffix = c(".dev", ".new")
  ) %>%
  # Step 2: Convert GPP means to integers
  mutate(
    mean_gpp_dev = as.integer(mean_gpp_dev),
    mean_gpp_ec  = as.integer(mean_gpp_ec),
    mean_gpp_qm  = as.integer(mean_gpp_qm)
  ) %>%
  # Step 3: Select and reorder columns
  select(
    koppen_class_name,
    n_sites_per_class,
    mean_gpp_ec,
    mean_gpp_dev,
    mean_gpp_qm,
    MAE.dev, RMSE.dev, NRMSE.dev,
    MAE.new, RMSE.new, NRMSE.new
  )

years_per_class <- rmse_dev_kopp_class_gpp$n_site_years_per_class

# Ensure n_site_years_per_class matches the total site-years
n_site_years_per_class = sum(rmse_new_kopp_class_gpp$n_site_years_per_class)

# Step 4: Prepare the overall summary row
overall_gpp_stats_for_kop_gramm2 <- overall_gpp_stats_for_kop_gramm %>%
  mutate(
    n_sites_per_class = sum(rmse_new_kopp_class_gpp$n_sites_per_class),
    koppen_class_name = "Overall"
  ) %>%
  # Align columns with rmse_total_by_kop_gpp
  select(names(rmse_total_by_kop_gpp))

# Step 5: Bind overall row with the class-level rows
rmse_total_by_kop_gpp2 <- bind_rows(
  overall_gpp_stats_for_kop_gramm2,
  rmse_total_by_kop_gpp
) %>%
  mutate(n_site_years_per_class = v <- append(rmse_dev_kopp_class_gpp$n_site_years_per_class, sum(rmse_dev_kopp_class_gpp$n_site_years_per_class), after = 0))


######################################################################################################
####################################### CREATE FIGURE B3 #############################################
######################################################################################################


#-----------------------------

# Define Functions

#-----------------------------

plot.dates.flux<- function(input_quincy_data, year_from, year_to, 
                           pheno_compare_full, input_flux) {
  gpp_out_plot <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = sos_quincy_out,
                  xmax = eos_quincy_out,
                  ymin = -Inf, ymax = Inf,
                  fill = "Output"),
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = sos_quincy_out,
                   color = "Output")) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = eos_quincy_out,
                   color = "Output")) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = sos_quincy_out,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(sos_quincy_out),
                   color = 'Output'),
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = eos_quincy_out,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(eos_quincy_out),
                   color = 'Output'),
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(name = 'Method',
                       breaks = "Output",
                       values = c('Output' = "black")) +
    scale_fill_manual(name = 'Method',
                      breaks = "Output",
                      values = c('Output' = "black")) +
    theme(legend.position = "top",
          axis.title.y = element_blank()) +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = TRS1.sos_quincy_elmore,
                  xmax = TRS1.eos_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "darkblue", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS1.sos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS1.eos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS1.sos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS1.sos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS1.eos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS1.eos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none") +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot_1 <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = TRS2.sos_quincy_elmore,
                  xmax = TRS2.eos_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "darkblue", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS2.sos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS2.eos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS2.sos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS2.sos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS2.eos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS2.eos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none") +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot_1a <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = TRS5.sos_quincy_elmore,
                  xmax = TRS5.eos_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "darkblue", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS5.sos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS5.eos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS5.sos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS5.sos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS5.eos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS5.eos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none") +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot_2 <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = DER.sos_quincy_elmore,
                  xmax = DER.eos_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "darkgreen", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = DER.sos_quincy_elmore),
               # linetype = 4,
               color = 'darkgreen') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = DER.eos_quincy_elmore),
               # linetype = 4,
               color = 'darkgreen') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = DER.sos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(DER.sos_quincy_elmore)),
               color = 'darkgreen',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = DER.eos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(DER.eos_quincy_elmore)),
               color = 'darkgreen',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # scale_x_date(date_labels = '%j') +
    theme_classic() +
    theme(legend.position = "none") +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  # ggtitle("QUINCY") +
  
  gpp_plot_3 <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = Greenup_quincy_elmore,
                  xmax = Dormancy_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "orange", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = Greenup_quincy_elmore),
               # linetype = 4,
               color = 'orange') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = Dormancy_quincy_elmore),
               # linetype = 4,
               color = 'orange') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = Greenup_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(Greenup_quincy_elmore)),
               color = 'orange',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = Dormancy_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(Dormancy_quincy_elmore)),
               color = 'orange',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none") +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot_4 <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = UD_quincy_elmore,
                  xmax = RD_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "#56B4E9", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = UD_quincy_elmore),
               # linetype = 4,
               color = '#56B4E9') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = RD_quincy_elmore),
               # linetype = 4,
               color = '#56B4E9') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = UD_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(UD_quincy_elmore)),
               color = '#56B4E9',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = RD_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(RD_quincy_elmore)),
               color = '#56B4E9',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none") +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  # fluxnet (with phen dates)
  gppflux_plot <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin =TRS1.sos_flux_elmore,
                  xmax = TRS1.eos_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "TRS1"), 
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS1.sos_flux_elmore,
                   color = 'TRS1')) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS1.eos_flux_elmore,
                   color = 'TRS1')) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS1.sos_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS1.sos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS1.eos_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS1.eos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.sos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.eos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
  theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1',
                                  'TRS2',
                                  'TRS5',
                                  'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  'TRS2' = "darkblue",
                                  'TRS5' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 
                                 'TRS2',
                                 'TRS5',
                                 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 'TRS2' = "darkblue",
                                 'TRS5' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    theme(legend.position = "top") +
    xlab("Time") +
    ylab("GPP (mol/day) - FLUX")
  
  gppflux_plot_1a <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin =TRS5.sos_flux_elmore,
                  xmax = TRS5.eos_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "TRS5"), 
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS5.sos_flux_elmore,
                   color = 'TRS5')) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS5.eos_flux_elmore,
                   color = 'TRS5')) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS5.sos_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS5.sos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS5.eos_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS5.eos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.sos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.eos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
  theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1',
                                  'TRS2',
                                  'TRS5',
                                  'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  'TRS2' = "darkblue",
                                  'TRS5' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 
                                 'TRS2',
                                 'TRS5',
                                 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 'TRS2' = "darkblue",
                                 'TRS5' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    theme(legend.position = "top") +
    xlab("Time") +
    ylab("GPP (mol/day) - FLUX")
  
  gppflux_plot_1 <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin =TRS2.sos_flux_elmore,
                  xmax = TRS2.eos_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "TRS2"), 
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS2.sos_flux_elmore,
                   color = 'TRS2')) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS2.eos_flux_elmore,
                   color = 'TRS2')) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS2.sos_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS2.sos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS2.eos_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS2.eos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.sos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.eos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
  theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1',
                                  'TRS2',
                                  'TRS5',
                                  'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  'TRS2' = "darkblue",
                                  'TRS5' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 
                                 'TRS2',
                                 'TRS5',
                                 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 'TRS2' = "darkblue",
                                 'TRS5' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    theme(legend.position = "top") +
    xlab("Time") +
    ylab("GPP (mol/day) - FLUX")
  
  gppflux_plot_2 <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = DER.sos_flux_elmore,
                  xmax = DER.eos_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "DER"), 
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = DER.sos_flux_elmore, 
                   color = "DER")) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = DER.eos_flux_elmore,
                   color = "DER")) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = DER.sos_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(DER.sos_flux_elmore)),
               color = 'darkgreen',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = DER.eos_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(DER.eos_flux_elmore)),
               color = 'darkgreen',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1', 'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    theme(legend.position = "top") +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    xlab("Time") +
    ylab("GPP (mol/day) - FLUX")
  
  gppflux_plot_3 <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = Greenup_flux_elmore,
                  xmax = Dormancy_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "Inflection"),
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = Greenup_flux_elmore,
                   color = 'Inflection')) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = Dormancy_flux_elmore,
                   color = 'Inflection')) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = Greenup_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(Greenup_flux_elmore)),
               color = 'orange',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = Dormancy_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(Dormancy_flux_elmore)),
               color = 'orange',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1', 'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    theme(legend.position = "top") +
    # ggtitle("Fluxnet") +
    xlab("Time") +
    ylab("GPP (mol/day) - FLUX")
  
  gppflux_plot_4 <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = UD_flux_elmore,
                  xmax = RD_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "Gu"),
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = UD_flux_elmore,
                   color = "Gu")) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = RD_flux_elmore,
                   color = "Gu")) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = UD_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(UD_flux_elmore)),
               color = '#56B4E9',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = RD_flux_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(RD_flux_elmore)),
               color = '#56B4E9',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1', 'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    theme(legend.position = "top") +
    xlab("Time") +
    ylab("GPP (mol/day) - FLUX")
  
  plot_1 <-
    grid.arrange(gppflux_plot,
                 gppflux_plot_1,
                 gppflux_plot_1a,
                 gppflux_plot_2,
                 gppflux_plot_3,
                 gppflux_plot_4,
                 gpp_plot,
                 gpp_plot_1,
                 gpp_plot_1a,
                 gpp_plot_2,
                 gpp_plot_3,
                 gpp_plot_4,
                 ncol = 6)
  
  compare_plot <-
    grid.arrange(gpp_out_plot,
                 plot_1,
                 widths = c(0.5, 2))
  
  return(compare_plot)
}





plot.dates.gcc <- function(input_quincy_data, year_from, year_to, 
                           pheno_compare_full, input_flux) {
  gpp_out_plot <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = sos_quincy_out,
                  xmax = eos_quincy_out,
                  ymin = -Inf, ymax = Inf,
                  fill = "Output"),
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = sos_quincy_out,
                   color = "Output")) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = eos_quincy_out,
                   color = "Output")) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = sos_quincy_out,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(sos_quincy_out),
                   color = 'Output'),
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = eos_quincy_out,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(eos_quincy_out),
                   color = 'Output'),
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(name = 'Method',
                       breaks = "Output",
                       values = c('Output' = "black")) +
    scale_fill_manual(name = 'Method',
                      breaks = "Output",
                      values = c('Output' = "black")) +
    theme(legend.position = "top",
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = TRS1.sos_quincy_elmore,
                  xmax = TRS1.eos_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "darkblue", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS1.sos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS1.eos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS1.sos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS1.sos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS1.eos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS1.eos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot_1 <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = TRS2.sos_quincy_elmore,
                  xmax = TRS2.eos_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "darkblue", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS2.sos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS2.eos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS2.sos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS2.sos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS2.eos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS2.eos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot_1a <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = TRS5.sos_quincy_elmore,
                  xmax = TRS5.eos_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "darkblue", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS5.sos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS5.eos_quincy_elmore),
               # linetype = 4,
               color = 'darkblue') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS5.sos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS5.sos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS5.eos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(TRS5.eos_quincy_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot_2 <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = DER.sos_quincy_elmore,
                  xmax = DER.eos_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "darkgreen", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = DER.sos_quincy_elmore),
               # linetype = 4,
               color = 'darkgreen') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = DER.eos_quincy_elmore),
               # linetype = 4,
               color = 'darkgreen') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = DER.sos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(DER.sos_quincy_elmore)),
               color = 'darkgreen',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = DER.eos_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(DER.eos_quincy_elmore)),
               color = 'darkgreen',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # scale_x_date(date_labels = '%j') +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  # ggtitle("QUINCY") +
  
  gpp_plot_3 <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = Greenup_quincy_elmore,
                  xmax = Dormancy_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "orange", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = Greenup_quincy_elmore),
               # linetype = 4,
               color = 'orange') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = Dormancy_quincy_elmore),
               # linetype = 4,
               color = 'orange') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = Greenup_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(Greenup_quincy_elmore)),
               color = 'orange',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = Dormancy_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(Dormancy_quincy_elmore)),
               color = 'orange',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  gpp_plot_4 <-
    input_quincy_data %>%
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = GPP,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = UD_quincy_elmore,
                  xmax = RD_quincy_elmore,
                  ymin = -Inf, ymax = Inf),
              fill = "#56B4E9", alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = UD_quincy_elmore),
               # linetype = 4,
               color = '#56B4E9') +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = RD_quincy_elmore),
               # linetype = 4,
               color = '#56B4E9') +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = UD_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(UD_quincy_elmore)),
               color = '#56B4E9',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = RD_quincy_elmore,
                   y = mean(input_quincy_data$GPP, na.rm = TRUE),
                   label = yday(RD_quincy_elmore)),
               color = '#56B4E9',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GPP (mol/ day) - QUINCY")
  
  # fluxnet (with phen dates)
  gppflux_plot <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = gcc_90,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin =TRS1.sos_flux_elmore,
                  xmax = TRS1.eos_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "TRS1"), 
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS1.sos_flux_elmore,
                   color = 'TRS1')) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS1.eos_flux_elmore,
                   color = 'TRS1')) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS1.sos_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(TRS1.sos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS1.eos_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(TRS1.eos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.sos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.eos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
  theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1',
                                  'TRS2',
                                  'TRS5',
                                  'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  'TRS2' = "darkblue",
                                  'TRS5' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 
                                 'TRS2',
                                 'TRS5',
                                 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 'TRS2' = "darkblue",
                                 'TRS5' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GCC (90th percentile)")
  
  gppflux_plot_1a <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = gcc_90,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin =TRS5.sos_flux_elmore,
                  xmax = TRS5.eos_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "TRS5"), 
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS5.sos_flux_elmore,
                   color = 'TRS5')) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS5.eos_flux_elmore,
                   color = 'TRS5')) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS5.sos_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(TRS5.sos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS5.eos_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(TRS5.eos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.sos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.eos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
  theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1',
                                  'TRS2',
                                  'TRS5',
                                  'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  'TRS2' = "darkblue",
                                  'TRS5' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 
                                 'TRS2',
                                 'TRS5',
                                 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 'TRS2' = "darkblue",
                                 'TRS5' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GCC (90th percentile)")
  
  gppflux_plot_1 <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = gcc_90,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin =TRS2.sos_flux_elmore,
                  xmax = TRS2.eos_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "TRS2"), 
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS2.sos_flux_elmore,
                   color = 'TRS2')) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = TRS2.eos_flux_elmore,
                   color = 'TRS2')) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS2.sos_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(TRS2.sos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = TRS2.eos_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(TRS2.eos_flux_elmore)),
               color = 'darkblue',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.sos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
    # geom_vline(data = pheno_compare_full %>%
    #              filter(year %in% 1995:1998),
    #            aes(xintercept = TRS2.eos_quincy_elmore),
    #            linetype = 4,
    #            color = 'blue') +
  theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1',
                                  'TRS2',
                                  'TRS5',
                                  'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  'TRS2' = "darkblue",
                                  'TRS5' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 
                                 'TRS2',
                                 'TRS5',
                                 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 'TRS2' = "darkblue",
                                 'TRS5' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GCC (90th percentile)")
  
  gppflux_plot_2 <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = gcc_90,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = DER.sos_flux_elmore,
                  xmax = DER.eos_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "DER"), 
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = DER.sos_flux_elmore, 
                   color = "DER")) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = DER.eos_flux_elmore,
                   color = "DER")) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = DER.sos_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(DER.sos_flux_elmore)),
               color = 'darkgreen',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = DER.eos_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(DER.eos_flux_elmore)),
               color = 'darkgreen',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1', 'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    xlab("Time") +
    ylab("GCC (90th percentile)")
  
  gppflux_plot_3 <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = gcc_90,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = Greenup_flux_elmore,
                  xmax = Dormancy_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "Inflection"),
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = Greenup_flux_elmore,
                   color = 'Inflection')) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = Dormancy_flux_elmore,
                   color = 'Inflection')) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = Greenup_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(Greenup_flux_elmore)),
               color = 'orange',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = Dormancy_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(Dormancy_flux_elmore)),
               color = 'orange',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    # linetype = 4) +
    theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1', 'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    # ggtitle("Fluxnet") +
    xlab("Time") +
    ylab("GCC (90th percentile)")
  
  gppflux_plot_4 <-
    input_flux %>% 
    filter(year %in% year_from:year_to) %>%
    ggplot() +
    geom_line(aes(x = date,
                  y = gcc_90,
                  group = year), color = "darkgray") +
    geom_rect(data = pheno_compare_full %>%
                filter(year %in% year_from:year_to),
              aes(xmin = UD_flux_elmore,
                  xmax = RD_flux_elmore,
                  ymin = -Inf, ymax = Inf, fill = "Gu"),
              alpha = 0.1) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = UD_flux_elmore,
                   color = "Gu")) +
    # linetype = 4) +
    geom_vline(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(xintercept = RD_flux_elmore,
                   color = "Gu")) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = UD_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(UD_flux_elmore)),
               color = '#56B4E9',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    geom_label(data = pheno_compare_full %>%
                 filter(year %in% year_from:year_to),
               aes(x = RD_flux_elmore,
                   y = mean(input_flux$gcc_90, na.rm = TRUE),
                   label = yday(RD_flux_elmore)),
               color = '#56B4E9',
               angle = 90,
               vjust = 1.5, size = 3,
               show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(name ='Method',
                       breaks = c('TRS1', 'DER',
                                  'Inflection',
                                  "Gu"),
                       values = c('TRS1' = "darkblue",
                                  "DER" = "darkgreen",
                                  'Inflection' = "orange",
                                  "Gu" = "#56B4E9")) +
    scale_fill_manual(name ='Method',
                      breaks = c('TRS1', 'DER',
                                 'Inflection',
                                 "Gu"),
                      values = c('TRS1' = "darkblue",
                                 "DER" = "darkgreen",
                                 'Inflection' = "orange",
                                 "Gu" = "#56B4E9")) +
    # scale_x_date(date_labels = '%j') + 
    # ggtitle("Fluxnet") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 25, hjust = 0.5)) +
    xlab("Time") +
    ylab("GCC (90th percentile)")
  
  plot_1 <-
    grid.arrange(gppflux_plot,
                 gppflux_plot_1,
                 gppflux_plot_1a,
                 gppflux_plot_2,
                 gppflux_plot_3,
                 gppflux_plot_4,
                 gpp_plot,
                 gpp_plot_1,
                 gpp_plot_1a,
                 gpp_plot_2,
                 gpp_plot_3,
                 gpp_plot_4,
                 ncol = 6)
  
  compare_plot <-
    grid.arrange(gpp_out_plot,
                 plot_1,
                 widths = c(0.5, 2))
  
  return(compare_plot)
}

# Save plots
save.plot <- 
  function(file_name, plot_name){
    ggplot2::ggsave(
      filename = file_name,
      plot = plot_name,
      width = 18,
      height = 10,
      dpi = 600)
  }





#------------------------------------------
# 2. Plot data ----
#------------------------------------------
(compare_cge_2_full_plot <-
   plot.dates.gcc(cge_quincy_c, 
                  2018,
                  2019, 
                  pheno_compare_cge,
                  cge_gcc))

(compare_cge_2_dev_full_plot <-
    plot.dates.gcc(cge_dev_quincy_c, 
                   2018,
                   2019, 
                   pheno_compare_cge_dev,
                   cge_gcc))

(compare_ustol_2_full_plot <-
    plot.dates.flux.2(ustol_quincy_c, 
                      2014,
                      2016, 
                      pheno_compare_ustol,
                      usics_flux_gpp))

(compare_ustol_dev_full_2_plot <-
    plot.dates.flux.2(ustol_dev_quincy_c, 
                      2014,
                      2016, 
                      pheno_compare_ustol_dev,
                      usics_flux_gpp))

(compare_esmag_2_full_plot <-
    plot.dates.flux.2(esmag_quincy_c, 
                      2004,
                      2006, 
                      pheno_compare_esmag,
                      esmag_flux_gpp))

(compare_esmag_dev_2_full_plot <-
    plot.dates.flux.2(esmag_dev_quincy_c, 
                      2004,
                      2006, 
                      pheno_compare_esmag_dev,
                      esmag_flux_gpp))

(compare_iedri_full_2_plot <-
    plot.dates.flux.2(iedri_quincy_c, 
                      2003,
                      2005, 
                      pheno_compare_iedri,
                      iedri_flux_gpp))

(compare_iedri_dev_full_2_plot <-
    plot.dates.flux.2(iedri_dev_quincy_c, 
                      2003,
                      2005, 
                      pheno_compare_iedri_dev,
                      iedri_flux_gpp))




# figures
figure_list <- 
  grep("_plot", names(.GlobalEnv), value = TRUE)

figure_list %>% 
  purrr::set_names() %>% 
  purrr::walk(
    .x = .,
    .f = ~ save.plot(
      file_name =  (paste0("path/to/output/file",.x,".png")),
      plot_name = get(.x)))






