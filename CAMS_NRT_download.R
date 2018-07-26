

library(RCurl)
library(stringr)
library(plyr)

setwd("/disk3/fkaragulian/ECMWF_forecasts/forecasts_downloads")
wd <- getwd()

date <- Sys.time()

url_list_data = "ftp://federico.karagulian:GI0cdFUW@dissemination.ecmwf.int/DATA/CAMS_NREALTIME/"
list_data = getURL(url_list_data, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE, ssl.verifypeer = FALSE) 
list_data = paste(url_list_data, strsplit(list_data, "\r*\n")[[1]], sep = "") 

# last element (date) of the list in order of time 
list_data <- sort(list_data)
n <- as.numeric(length(list_data))
list_data[n]
folder <- unlist(list_data[n])

# isolate folder name
folder <- str_sub(folder, start = 80, end = -1)

url = paste("ftp://federico.karagulian:GI0cdFUW@dissemination.ecmwf.int/DATA/CAMS_NREALTIME/",folder, "/", sep = "")

filenames = getURL(url, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE, ssl.verifypeer = FALSE) 
filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "") 

# create a new working directory for each download date
dir.create(folder)

# PM2.5
setwd(paste0(wd,"/",folder))
dir.create("PM25")
setwd(paste0(wd,"/",folder,"/PM25"))
filenames_pm25 <- unlist(str_extract_all(filenames, ".+(_pm2p5.nc$)"))
filenames_pm25 <- sort(filenames_pm25)
filenames_pm25 <- filenames_pm25[109:120]
# start downloading data in the main directory -----------------------

# mapply(download.file, filenames_pm25,basename(filenames_pm25), method = 'curl') 
mapply(download.file, filenames_pm25,basename(filenames_pm25)) 

# PM10
setwd(paste0(wd,"/",folder))
dir.create("PM10")
setwd(paste0(wd,"/",folder,"/PM10"))
filenames_pm10 <- unlist(str_extract_all(filenames, ".+(_pm10.nc$)"))
filenames_pm10 <- sort(filenames_pm10)
filenames_pm10 <- filenames_pm10[109:120]
mapply(download.file, filenames_pm10,basename(filenames_pm10)) 

# NO2
setwd(paste0(wd,"/",folder))
dir.create("NO2")
setwd(paste0(wd,"/",folder,"/NO2"))
filenames_no2 <- unlist(str_extract_all(filenames, ".+(_tcno2.nc$)"))
filenames_no2 <- sort(filenames_no2)
filenames_no2 <- filenames_no2[110:121]
mapply(download.file, filenames_no2,basename(filenames_no2)) 

# SO2
setwd(paste0(wd,"/",folder))
dir.create("SO2")
setwd(paste0(wd,"/",folder,"/SO2"))
filenames_so2 <- unlist(str_extract_all(filenames, ".+(_tcso2.nc$)"))
filenames_so2 <- sort(filenames_so2)
filenames_so2 <- filenames_so2[110:121]
mapply(download.file, filenames_so2,basename(filenames_so2)) 


# CO
setwd(paste0(wd,"/",folder))
dir.create("CO")
setwd(paste0(wd,"/",folder,"/CO"))
filenames_co <- unlist(str_extract_all(filenames, ".+(_tcco.nc$)"))
filenames_co <- sort(filenames_co)
filenames_co <- filenames_co[110:121]
mapply(download.file, filenames_co,basename(filenames_co)) 


# AOD 550 nm
setwd(paste0(wd,"/",folder))
dir.create("AOD")
setwd(paste0(wd,"/",folder,"/AOD"))
filenames_aod <- unlist(str_extract_all(filenames, ".+(_aod550.nc$)"))
filenames_aod <- sort(filenames_aod)
filenames_aod <- filenames_aod[109:120]
mapply(download.file, filenames_aod,basename(filenames_aod)) 


# sea salt AOD 550 nm
setwd(paste0(wd,"/",folder))
dir.create("SS_AOD")
setwd(paste0(wd,"/",folder,"/SS_AOD"))
filenames_ssaod <- unlist(str_extract_all(filenames, ".+(_ssaod550.nc$)"))
filenames_ssaod <- sort(filenames_ssaod)
filenames_ssaod <- filenames_ssaod[109:120]
mapply(download.file, filenames_ssaod,basename(filenames_ssaod)) 



# dust AOD 550 nm
setwd(paste0(wd,"/",folder))
dir.create("DU_AOD")
setwd(paste0(wd,"/",folder,"/DU_AOD"))
filenames_duaod <- unlist(str_extract_all(filenames, ".+(_duaod550.nc$)"))
filenames_duaod <- sort(filenames_duaod)
filenames_duaod <- filenames_duaod[109:120]
mapply(download.file, filenames_duaod,basename(filenames_duaod)) 

# Organic Matter AOD 550 nm
setwd(paste0(wd,"/",folder))
dir.create("OM_AOD")
setwd(paste0(wd,"/",folder,"/OM_AOD"))
filenames_omaod <- unlist(str_extract_all(filenames, ".+(_omaod550.nc$)"))
filenames_omaod <- sort(filenames_omaod)
filenames_omaod <- filenames_omaod[109:120]
mapply(download.file, filenames_omaod,basename(filenames_omaod)) 

# Black carbon AOD 550 nm
setwd(paste0(wd,"/",folder))
dir.create("BC_AOD")
setwd(paste0(wd,"/",folder,"/BC_AOD"))
filenames_bcaod <- unlist(str_extract_all(filenames, ".+(_bcaod550.nc$)"))
filenames_bcaod <- sort(filenames_bcaod)
filenames_bcaod <- filenames_bcaod[109:120]
mapply(download.file, filenames_bcaod,basename(filenames_bcaod)) 

# HNO3 550 nm
setwd(paste0(wd,"/",folder))
dir.create("HNO3")
setwd(paste0(wd,"/",folder,"/HNO3"))
filenames_hno3 <- unlist(str_extract_all(filenames, ".+(_tc_hno3.nc$)"))
filenames_hno3 <- sort(filenames_hno3)
filenames_hno3 <- filenames_hno3[110:121]
mapply(download.file, filenames_hno3,basename(filenames_hno3)) 


# NO
setwd(paste0(wd,"/",folder))
dir.create("NO")
setwd(paste0(wd,"/",folder,"/NO"))
filenames_no <- unlist(str_extract_all(filenames, ".+(_tc_no.nc$)"))
filenames_no <- sort(filenames_no)
filenames_no <- filenames_no[110:121]
mapply(download.file, filenames_no,basename(filenames_no)) 

