
library(readr)
library(dplyr)

# bind all regridded satellite data from ECMWF

# load ECMWF data
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/ECMWF_MODIS_csv_40km")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2015/ECMWF_MODIS_csv_40km")
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2014/ECMWF_MODIS_csv_40km")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2013/ECMWF_MODIS_csv_40km")


filenames <- list.files(pattern = "\\.csv$")

### Function to bind all AOD ########################

DATE = NULL
LAT = NULL
LON = NULL
AOD_MODIS = NULL
AOD_ECMWF = NULL
PM25_ECMWF = NULL
PM10_ECMWF = NULL
DUST_ECMWF = NULL
SALT_ECMWF = NULL
BC_ECMWF = NULL
OM_ECMWF = NULL
SO4_ECMWF = NULL



# filenames <- filenames[3]

  for (i in 1:length(filenames)) {

    # date <- read_csv(filenames[i])[,1]
    # lon <- read_csv(filenames[i])[,2]
    # lat <- read_csv(filenames[i])[,3]
    # aod_MODIS <- read_csv(filenames[i])[,4]
    # aod_ECMWF <- read_csv(filenames[i])[,5]
    # pm25_ECMWF <- read_csv(filenames[i])[,6]
    # pm10_ECMWF <- read_csv(filenames[i])[,7]
    # dust_ECMWF <- read_csv(filenames[i])[,8]
    # salt_ECMWF <- read_csv(filenames[i])[,9]
    # bc_ECMWF <- read_csv(filenames[i])[,10]
    # om_ECMWF <- read_csv(filenames[i])[,11]
    # so4_ECMWF <- read_csv(filenames[i])[,12]
    
    
    date <- read_csv(filenames[i])[,1]
    lon <- read_csv(filenames[i])[,2]
    lat <- read_csv(filenames[i])[,3]
    aod_MODIS <- read_csv(filenames[i])[,4]
    aod_ECMWF <- read_csv(filenames[i])[,5]
    dust_ECMWF <- read_csv(filenames[i])[,6]
    salt_ECMWF <- read_csv(filenames[i])[,7]
    bc_ECMWF <- read_csv(filenames[i])[,8]
    om_ECMWF <- read_csv(filenames[i])[,9]
    so4_ECMWF <- read_csv(filenames[i])[,10]
    
    DATE = rbind(DATE, data.frame(date))
    LON = rbind(LON, data.frame(lon))
    LAT = rbind(LAT, data.frame(lat))
    AOD_MODIS = rbind(AOD_MODIS, data.frame(aod_MODIS))
    AOD_ECMWF = rbind(AOD_ECMWF, data.frame(aod_ECMWF))
    # PM25_ECMWF = rbind(PM25_ECMWF, data.frame(pm25_ECMWF))
    # PM10_ECMWF = rbind(PM10_ECMWF, data.frame(pm10_ECMWF))
    DUST_ECMWF = rbind(DUST_ECMWF, data.frame(dust_ECMWF))
    SALT_ECMWF = rbind(SALT_ECMWF, data.frame(salt_ECMWF))
    BC_ECMWF = rbind(BC_ECMWF, data.frame(bc_ECMWF))
    OM_ECMWF = rbind(OM_ECMWF, data.frame(om_ECMWF))
    SO4_ECMWF = rbind(SO4_ECMWF, data.frame(so4_ECMWF))
    
}
  

AOD_data <- cbind(DATE, LON, LAT, 
                  AOD_MODIS,
                  AOD_ECMWF,
                  # PM25_ECMWF,
                  # PM10_ECMWF,
                  DUST_ECMWF,
                  SALT_ECMWF,
                  BC_ECMWF,
                  OM_ECMWF,
                  SO4_ECMWF)

# write_csv(AOD_data, "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/AOD_2016_MODIS_all_40km.csv")
# write_csv(AOD_data, "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2015/AOD_2015_MODIS_all_40km.csv")
write_csv(AOD_data, "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2014/AOD_2014_MODIS_all_40km.csv")
# write_csv(AOD_data, "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2013/AOD_2013_MODIS_all_40km.csv")


