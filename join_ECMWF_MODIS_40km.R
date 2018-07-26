
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(splines)
library(plyr)
library(stringr)

# ECMWF data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2013/all_ECMWF")
 setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2014/all_ECMWF")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2015/all_ECMWF")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/all_ECMWF")

# create a new directory for ECMWF & MODIS joint files
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2013/ECMWF_MODIS_csv_40km")
 dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2014/ECMWF_MODIS_csv_40km")
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2015/ECMWF_MODIS_csv_40km")
#  dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/ECMWF_MODIS_csv_40km")

 
# MODIS data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
# new_dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2013/ECMWF_MODIS_csv_40km"
 new_dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2014/ECMWF_MODIS_csv_40km"
# new_dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2015/ECMWF_MODIS_csv_40km"
# new_dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/ECMWF_MODIS_csv_40km"

filenames_ECMWF <- list.files(pattern = "\\.csv$")

# MODIS processed data at 40km
# dir_MODIS_AOD <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/csv_40km"
 dir_MODIS_AOD <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2014_MODIS_processed/csv_40km"
# dir_MODIS_AOD <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/csv_40km"
# dir_MODIS_AOD <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed/csv_40km"

filenames_MODIS_AOD <- list.files(path = dir_MODIS_AOD, pattern = "\\.csv$")

# i <- 1
# filenames_ECMWF <- filenames_ECMWF[i]

for (i in 1: length(filenames_ECMWF)) {

ECMWF <- read_csv(filenames_ECMWF[i])
MODIS_AOD <- read_csv(paste0(dir_MODIS_AOD, "/", filenames_MODIS_AOD[i]))
date <- str_sub(filenames_ECMWF[i], start = 1, end = -5)

# joint data on with same Lon and Lat

ECMWF_MODIS_joint <- MODIS_AOD %>%
  merge(ECMWF, by = c("Date", "Lon", "Lat"))
# colnames(ECMWF_MODIS_joint) <- c("Date", "Lon", "Lat", "AOD_MODIS", "AOD_ECMWF",
#                                  "PM25_ECMWF", "PM10_ECMWF", "DUST_ECMWF",
#                                   "SALT_ECMWF", "BC_ECMWF",
#                                  "OM_ECMWF", "SO4_ECMWF")


colnames(ECMWF_MODIS_joint) <- c("Date", "Lon", "Lat", "AOD_MODIS", "AOD_ECMWF",
                                 "DUST_ECMWF",
                                 "SALT_ECMWF", "BC_ECMWF",
                                 "OM_ECMWF", "SO4_ECMWF")

write_csv(ECMWF_MODIS_joint, paste0(path = new_dir,"/", date,".csv"))

}




