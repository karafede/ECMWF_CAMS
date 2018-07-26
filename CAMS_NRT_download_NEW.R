

library(RCurl)
library(stringr)
library(plyr)
library(RNetCDF)
library(rgdal)
library(ncdf4)
library(raster)
library(stringr)
library(lubridate)
library(gstat)
library(ncdf.tools)


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
# setwd(paste0(wd,"/",folder))
# dir.create("PM10")
# setwd(paste0(wd,"/",folder,"/PM10"))
# filenames_pm10 <- unlist(str_extract_all(filenames, ".+(_pm10.nc$)"))
# filenames_pm10 <- sort(filenames_pm10)
# filenames_pm10 <- filenames_pm10[109:120]
# mapply(download.file, filenames_pm10,basename(filenames_pm10)) 

# NO2
# setwd(paste0(wd,"/",folder))
# dir.create("NO2")
# setwd(paste0(wd,"/",folder,"/NO2"))
# filenames_no2 <- unlist(str_extract_all(filenames, ".+(_tcno2.nc$)"))
# filenames_no2 <- sort(filenames_no2)
# filenames_no2 <- filenames_no2[110:121]
# mapply(download.file, filenames_no2,basename(filenames_no2)) 

# SO2
# setwd(paste0(wd,"/",folder))
# dir.create("SO2")
# setwd(paste0(wd,"/",folder,"/SO2"))
# filenames_so2 <- unlist(str_extract_all(filenames, ".+(_tcso2.nc$)"))
# filenames_so2 <- sort(filenames_so2)
# filenames_so2 <- filenames_so2[110:121]
# mapply(download.file, filenames_so2,basename(filenames_so2)) 


# CO
# setwd(paste0(wd,"/",folder))
# dir.create("CO")
# setwd(paste0(wd,"/",folder,"/CO"))
# filenames_co <- unlist(str_extract_all(filenames, ".+(_tcco.nc$)"))
# filenames_co <- sort(filenames_co)
# filenames_co <- filenames_co[110:121]
# mapply(download.file, filenames_co,basename(filenames_co)) 

# O3
# setwd(paste0(wd,"/",folder))
# dir.create("O3")
# setwd(paste0(wd,"/",folder,"/O3"))
# filenames_O3 <- unlist(str_extract_all(filenames, ".+(_gtco3.nc$)"))
# filenames_O3 <- sort(filenames_O3)
# filenames_O3 <- filenames_O3[110:121]
# mapply(download.file, filenames_O3,basename(filenames_O3)) 


# # AOD 550 nm
# setwd(paste0(wd,"/",folder))
# dir.create("AOD")
# setwd(paste0(wd,"/",folder,"/AOD"))
# filenames_aod <- unlist(str_extract_all(filenames, ".+(_aod550.nc$)"))
# filenames_aod <- sort(filenames_aod)
# filenames_aod <- filenames_aod[109:120]
# mapply(download.file, filenames_aod,basename(filenames_aod)) 
# 
# 
# # sea salt AOD 550 nm
# setwd(paste0(wd,"/",folder))
# dir.create("SS_AOD")
# setwd(paste0(wd,"/",folder,"/SS_AOD"))
# filenames_ssaod <- unlist(str_extract_all(filenames, ".+(_ssaod550.nc$)"))
# filenames_ssaod <- sort(filenames_ssaod)
# filenames_ssaod <- filenames_ssaod[109:120]
# mapply(download.file, filenames_ssaod,basename(filenames_ssaod)) 
# 
# 
# 
# # dust AOD 550 nm
# setwd(paste0(wd,"/",folder))
# dir.create("DU_AOD")
# setwd(paste0(wd,"/",folder,"/DU_AOD"))
# filenames_duaod <- unlist(str_extract_all(filenames, ".+(_duaod550.nc$)"))
# filenames_duaod <- sort(filenames_duaod)
# filenames_duaod <- filenames_duaod[109:120]
# mapply(download.file, filenames_duaod,basename(filenames_duaod)) 
# 
# # Organic Matter AOD 550 nm
# setwd(paste0(wd,"/",folder))
# dir.create("OM_AOD")
# setwd(paste0(wd,"/",folder,"/OM_AOD"))
# filenames_omaod <- unlist(str_extract_all(filenames, ".+(_omaod550.nc$)"))
# filenames_omaod <- sort(filenames_omaod)
# filenames_omaod <- filenames_omaod[109:120]
# mapply(download.file, filenames_omaod,basename(filenames_omaod)) 
# 
# Black carbon AOD 550 nm
# setwd(paste0(wd,"/",folder))
# dir.create("BC_AOD")
# setwd(paste0(wd,"/",folder,"/BC_AOD"))
# filenames_bcaod <- unlist(str_extract_all(filenames, ".+(_bcaod550.nc$)"))
# filenames_bcaod <- sort(filenames_bcaod)
# filenames_bcaod <- filenames_bcaod[109:120]
# mapply(download.file, filenames_bcaod,basename(filenames_bcaod))
# 
# # HNO3 550 nm
# setwd(paste0(wd,"/",folder))
# dir.create("HNO3")
# setwd(paste0(wd,"/",folder,"/HNO3"))
# filenames_hno3 <- unlist(str_extract_all(filenames, ".+(_tc_hno3.nc$)"))
# filenames_hno3 <- sort(filenames_hno3)
# filenames_hno3 <- filenames_hno3[110:121]
# mapply(download.file, filenames_hno3,basename(filenames_hno3)) 
# 
# 
# # NO
# setwd(paste0(wd,"/",folder))
# dir.create("NO")
# setwd(paste0(wd,"/",folder,"/NO"))
# filenames_no <- unlist(str_extract_all(filenames, ".+(_tc_no.nc$)"))
# filenames_no <- sort(filenames_no)
# filenames_no <- filenames_no[110:121]
# mapply(download.file, filenames_no,basename(filenames_no)) 


#####################################################################################
## get only PM2.5 ###################################################################

# setwd("D:/ECMWF_data_download")
setwd("/disk3/fkaragulian/ECMWF_forecasts/forecasts_downloads/")


# dir <- "D:/MODIS_AOD/UAE_boundary"
dir <- "/disk3/fkaragulian/MODIS_AOD/UAE_boundary"
### shapefile for UAE
shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")

# dir <- "D:/Dust_Event_UAE_2015/WRFChem_domain"

# larger WRF domain
shp_WRF <- readOGR(dsn = dir, layer = "domain_d02_4km_WRFChem_small")
# plot(shp_WRF)
# plot(shp_UAE, add=TRUE, lwd=1)



folders <- list.files()
# get the most recent directory
folders <- sort(folders, decreasing = T)

wd <- getwd()
# go to the PM25 folder
setwd(paste0(wd,"/",folders[1],"/PM25"))


patt<- ".nc"
filenames <- list.files(pattern = patt)
# filenames <- filenames[1]


# gerate a time sequence of 12 hours
TS <- 1:12

# i <- 1

for(i in 1:length(filenames)) {
  ECMWF_file <- open.nc(filenames[i])
  ECMWF_file <- read.nc(ECMWF_file)
  name_vari <- names(ECMWF_file)
  
  
  
  datetime <- ECMWF_file[3] 
  datetime <- as.numeric(datetime)
  datetime <- convertDateNcdf2R(datetime,units="hours",origin = as.POSIXct("1900-01-01 00:00:01", tz = "UTC"), 
                                time.format = c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S", 
                                                "%Y-%m-%d %H:%M", "%Y-%m-%d %Z %H:%M", "%Y-%m-%d %Z %H:%M:%S"))
  datetime <- as.character(datetime)
  time <- as.numeric(str_sub(filenames[i], start = 23, end = -30))
  # time <- 0
  # time <- 12
  if (time == 0) {
    time <- 11
    
    year <- str_sub(datetime, start = 0, end = -16)
    month <- str_sub(datetime, start = 6, end = -13)
    day <- str_sub(datetime, start = 9, end = -10) 
    
    date <- paste0(year,"_",month,"_", day)
    
  } else {
    time <- 0

    # year <- str_sub(datetime, start = 0, end = -7)
    # month <- str_sub(datetime, start = 6, end = -4)
    # day <- str_sub(datetime, start = 9, end = -1) 
    
    year <- str_sub(datetime, start = 0, end = -16)
    month <- str_sub(datetime, start = 6, end = -13)
    day <- str_sub(datetime, start = 9, end = -10) 
    
    date <- paste0(year,"_",month,"_", day)
    
  }
  
  
  PM25 <- (ECMWF_file[4])    #  only one variable (pm2p5) == var = 4
  names(PM25)<- "xxyyzz"
  PM25<- (PM25$xxyyzz)
  LON <- (ECMWF_file$longitude) -180
  # LON_West <- LON[1:451] -180
  # LON_East <- LON[452:900] -180
  # LON <- c(LON_West, LON_East)
  LAT <- ECMWF_file$latitude
  
  
  xmn= min(LON)
  xmx=max(LON)
  ymn=min(LAT)
  ymx=max(LAT)

  MMM <-  PM25[ , ]*1e+09
  MMM <- MMM[ c((451:1),(900:452)), ]
  MMM <- MMM[nrow(MMM):1, ]
  MMM <-  t(MMM[ , ])  # map is upside down 
  r <- raster(MMM, xmn, xmx, ymn,  ymx, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
 # plot(r)
  name_time <- time + TS[i]
  name_time <- str_pad(name_time, 4, pad="0")
  name_time <- str_sub(name_time, start = 3, end = -1)
  name_time <- paste0(name_time, "_00")
  
  # crop over the UAE area
   r <- crop(r, extent(shp_WRF))
   r <- mask(r, shp_WRF)
#   plot(r)
   
   # kirging
   ### Extract points from raster tiff ############################################
   
   r_pts <- rasterToPoints(r)
   colnames(r_pts) <- c("lon", "lat", "PM25")
   r_pts <- as.data.frame(r_pts) 
   
   # converta data into a spatialpolygon dataframe------------------------------------
   
   r_pts$x <- r_pts$lon
   r_pts$y <- r_pts$lat
   
   coordinates(r_pts) = ~x + y  ## Set spatial coordinates to create a Spatial object:
   
   
   # make a variogram----------------------------------------------------------------
   
   vargram_ECMWF <- variogram(PM25 ~ 1, r_pts) # calculates sample variogram values 
   
   # fite the variogram
   vargram_ECMWF_fit  <- fit.variogram(vargram_ECMWF, fit.ranges = FALSE, fit.sills = FALSE,
                                     vgm(1, c("Sph", "Exp", "Mat")), fit.kappa = TRUE)
   
   # plot(vargram_ECMWF, vargram_ECMWF_fit) # plot the sample values, along with the fit model
   
   e <- extent(shp_WRF)
   
   # make a regular empty grid
   x.range <- as.numeric(c(e@xmin, e@xmax))  # min/max longitude of the interpolation area
   y.range <- as.numeric(c(e@ymin, e@ymax))  # min/max latitude of the interpolation area
   
   ## grid at 10km resolution
   grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1),
                      y = seq(from = y.range[1], to = y.range[2], by = 0.1))  # expand points to grid
   coordinates(grd) <- ~x + y
   gridded(grd) <- TRUE
   
   f.1 <- as.formula(Precip_in ~ X + Y)
   # perform kriging
   dat.krg <- gstat::krige(PM25 ~ 1, r_pts, grd, vargram_ECMWF_fit, nmax = 50)
   r <- raster(dat.krg)
   projection(r) <- CRS("+proj=longlat +datum=WGS84")
   # r <- crop(r, extent(shp_UAE))
   # r <- mask(r, shp_UAE)
   
   # plot(r)
   # plot(shp_UAE, add = T, lwd =1)
   
   
   
   #########-------------------------------------------------------------------------
   # load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)
   
   PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   PM25_2015_GRW <- raster("D:/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
   projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")
   
   ### crop raster over the UAE shp file  ###############################
   PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
   PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)
   
   ###################
   
   ## make resolution of MODIS-data as the one of GWR-------------------------------
   ECMWF_tif_1km = projectRaster(r, PM25_2015_GRW)
   ECMWF_tif_1km <- crop(ECMWF_tif_1km, extent(shp_UAE))
   ECMWF_tif_1km <- mask(ECMWF_tif_1km, shp_UAE)
   
 #  plot(ECMWF_tif_1km)
   
#  writeRaster(ECMWF_tif_1km, paste0("D:/ECMWF_data_download/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
  writeRaster(ECMWF_tif_1km, paste0(wd,"/",folders[1],"/PM25/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)

}

#####################################################################################
#####################################################################################
#####################################################################################
## get SO2 ##########################################################################

# # setwd("D:/ECMWF_data_download")
# setwd("/disk3/fkaragulian/ECMWF_forecasts/forecasts_downloads/")
# 
# 
# # dir <- "D:/MODIS_AOD/UAE_boundary"
# dir <- "/disk3/fkaragulian/MODIS_AOD/UAE_boundary"
# ### shapefile for UAE
# shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")
# 
# # dir <- "D:/Dust_Event_UAE_2015/WRFChem_domain"
# 
# # larger WRF domain
# shp_WRF <- readOGR(dsn = dir, layer = "domain_d02_4km_WRFChem_small")
# # plot(shp_WRF)
# # plot(shp_UAE, add=TRUE, lwd=1)
# 
# 
# 
# folders <- list.files()
# # get the most recent directory
# folders <- sort(folders, decreasing = T)
# 
# wd <- getwd()
# # go to the PM25 folder
# setwd(paste0(wd,"/",folders[1],"/SO2"))
# 
# 
# patt<- ".nc"
# filenames <- list.files(pattern = patt)
# # filenames <- filenames[1]
# 
# 
# # gerate a time sequence of 12 hours
# TS <- 1:12
# 
# # i <- 1
# 
# for(i in 1:length(filenames)) {
#   ECMWF_file <- open.nc(filenames[i])
#   ECMWF_file <- read.nc(ECMWF_file)
#   name_vari <- names(ECMWF_file)
#   
#   
#   
#   datetime <- ECMWF_file[3] 
#   datetime <- as.numeric(datetime)
#   datetime <- convertDateNcdf2R(datetime,units="hours",origin = as.POSIXct("1900-01-01 00:00:01", tz = "UTC"), 
#                                 time.format = c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S", 
#                                                 "%Y-%m-%d %H:%M", "%Y-%m-%d %Z %H:%M", "%Y-%m-%d %Z %H:%M:%S"))
#   datetime <- as.character(datetime)
#   time <- as.numeric(str_sub(filenames[i], start = 23, end = -30))
#   # time <- 0
#   # time <- 12
#   if (time == 0) {
#     time <- 11
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   } else {
#     time <- 0
#     
#     # year <- str_sub(datetime, start = 0, end = -7)
#     # month <- str_sub(datetime, start = 6, end = -4)
#     # day <- str_sub(datetime, start = 9, end = -1) 
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   }
#   
#   
#   SO2 <- (ECMWF_file[4])    #  only one variable (tcso2) == var = 4
#   names(SO2)<- "xxyyzz"
#   SO2 <- (SO2$xxyyzz)
#   LON <- (ECMWF_file$longitude) -180
#   # LON_West <- LON[1:451] -180
#   # LON_East <- LON[452:900] -180
#   # LON <- c(LON_West, LON_East)
#   LAT <- ECMWF_file$latitude
#   
#   
#   xmn= min(LON)
#   xmx=max(LON)
#   ymn=min(LAT)
#   ymx=max(LAT)
#   
#   MMM <-  SO2[ , ]*1e+09
#   MMM <- MMM[ c((451:1),(900:452)), ]
#   MMM <- MMM[nrow(MMM):1, ]
#   MMM <-  t(MMM[ , ])  # map is upside down 
#   r <- raster(MMM, xmn, xmx, ymn,  ymx, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#   # plot(r)
#   name_time <- time + TS[i]
#   name_time <- str_pad(name_time, 4, pad="0")
#   name_time <- str_sub(name_time, start = 3, end = -1)
#   name_time <- paste0(name_time, "_00")
#   
#   # crop over the UAE area
#   r <- crop(r, extent(shp_WRF))
#   r <- mask(r, shp_WRF)
#   #   plot(r)
#   
#   # kirging
#   ### Extract points from raster tiff ############################################
#   
#   r_pts <- rasterToPoints(r)
#   colnames(r_pts) <- c("lon", "lat", "SO2")
#   r_pts <- as.data.frame(r_pts) 
#   
#   # converta data into a spatialpolygon dataframe------------------------------------
#   
#   r_pts$x <- r_pts$lon
#   r_pts$y <- r_pts$lat
#   
#   coordinates(r_pts) = ~x + y  ## Set spatial coordinates to create a Spatial object:
#   
#   
#   # make a variogram----------------------------------------------------------------
#   
#   vargram_ECMWF <- variogram(SO2 ~ 1, r_pts) # calculates sample variogram values 
#   
#   # fite the variogram
#   vargram_ECMWF_fit  <- fit.variogram(vargram_ECMWF, fit.ranges = FALSE, fit.sills = FALSE,
#                                       vgm(1, c("Sph", "Exp", "Mat")), fit.kappa = TRUE)
#   
#   # plot(vargram_ECMWF, vargram_ECMWF_fit) # plot the sample values, along with the fit model
#   
#   e <- extent(shp_WRF)
#   
#   # make a regular empty grid
#   x.range <- as.numeric(c(e@xmin, e@xmax))  # min/max longitude of the interpolation area
#   y.range <- as.numeric(c(e@ymin, e@ymax))  # min/max latitude of the interpolation area
#   
#   ## grid at 10km resolution
#   grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1),
#                      y = seq(from = y.range[1], to = y.range[2], by = 0.1))  # expand points to grid
#   coordinates(grd) <- ~x + y
#   gridded(grd) <- TRUE
#   
#   f.1 <- as.formula(Precip_in ~ X + Y)
#   # perform kriging
#   dat.krg <- gstat::krige(SO2 ~ 1, r_pts, grd, vargram_ECMWF_fit, nmax = 50)
#   r <- raster(dat.krg)
#   projection(r) <- CRS("+proj=longlat +datum=WGS84")
#   # r <- crop(r, extent(shp_UAE))
#   # r <- mask(r, shp_UAE)
#   
#   # plot(r)
#   # plot(shp_UAE, add = T, lwd =1)
#   
#   
#   
#   #########-------------------------------------------------------------------------
#   # load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)
#   
#   PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   #   PM25_2015_GRW <- raster("D:/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")
#   
#   ### crop raster over the UAE shp file  ###############################
#   PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
#   PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)
#   
#   ###################
#   
#   ## make resolution of MODIS-data as the one of GWR-------------------------------
#   ECMWF_tif_1km = projectRaster(r, PM25_2015_GRW)
#   ECMWF_tif_1km <- crop(ECMWF_tif_1km, extent(shp_UAE))
#   ECMWF_tif_1km <- mask(ECMWF_tif_1km, shp_UAE)
#   
#   #  plot(ECMWF_tif_1km)
#   
#   #  writeRaster(ECMWF_tif_1km, paste0("D:/ECMWF_data_download/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   writeRaster(ECMWF_tif_1km, paste0(wd,"/",folders[1],"/SO2/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   
# }
# 
# 
# #####################################################################################
# #####################################################################################
# ## get O3 ###########################################################################
# #####################################################################################
# #####################################################################################
# 
# # setwd("D:/ECMWF_data_download")
# setwd("/disk3/fkaragulian/ECMWF_forecasts/forecasts_downloads/")
# 
# 
# # dir <- "D:/MODIS_AOD/UAE_boundary"
# dir <- "/disk3/fkaragulian/MODIS_AOD/UAE_boundary"
# ### shapefile for UAE
# shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")
# 
# # dir <- "D:/Dust_Event_UAE_2015/WRFChem_domain"
# 
# # larger WRF domain
# shp_WRF <- readOGR(dsn = dir, layer = "domain_d02_4km_WRFChem_small")
# # plot(shp_WRF)
# # plot(shp_UAE, add=TRUE, lwd=1)
# 
# 
# 
# folders <- list.files()
# # get the most recent directory
# folders <- sort(folders, decreasing = T)
# 
# wd <- getwd()
# # go to the PM25 folder
# setwd(paste0(wd,"/",folders[1],"/O3"))
# 
# 
# patt<- ".nc"
# filenames <- list.files(pattern = patt)
# # filenames <- filenames[1]
# 
# 
# # gerate a time sequence of 12 hours
# TS <- 1:12
# 
# # i <- 1
# 
# for(i in 1:length(filenames)) {
#   ECMWF_file <- open.nc(filenames[i])
#   ECMWF_file <- read.nc(ECMWF_file)
#   name_vari <- names(ECMWF_file)
#   
#   
#   
#   datetime <- ECMWF_file[3] 
#   datetime <- as.numeric(datetime)
#   datetime <- convertDateNcdf2R(datetime,units="hours",origin = as.POSIXct("1900-01-01 00:00:01", tz = "UTC"), 
#                                 time.format = c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S", 
#                                                 "%Y-%m-%d %H:%M", "%Y-%m-%d %Z %H:%M", "%Y-%m-%d %Z %H:%M:%S"))
#   datetime <- as.character(datetime)
#   time <- as.numeric(str_sub(filenames[i], start = 23, end = -30))
#   # time <- 0
#   # time <- 12
#   if (time == 0) {
#     time <- 11
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   } else {
#     time <- 0
#     
#     # year <- str_sub(datetime, start = 0, end = -7)
#     # month <- str_sub(datetime, start = 6, end = -4)
#     # day <- str_sub(datetime, start = 9, end = -1) 
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   }
#   
#   
#   O3 <- (ECMWF_file[4])    #  only one variable (gtco3) == var = 4
#   names(O3)<- "xxyyzz"
#   O3 <- (O3$xxyyzz)
#   LON <- (ECMWF_file$longitude) -180
#   # LON_West <- LON[1:451] -180
#   # LON_East <- LON[452:900] -180
#   # LON <- c(LON_West, LON_East)
#   LAT <- ECMWF_file$latitude
#   
#   
#   xmn= min(LON)
#   xmx=max(LON)
#   ymn=min(LAT)
#   ymx=max(LAT)
#   
#   MMM <-  O3[ , ]*1e+09
#   MMM <- MMM[ c((451:1),(900:452)), ]
#   MMM <- MMM[nrow(MMM):1, ]
#   MMM <-  t(MMM[ , ])  # map is upside down 
#   r <- raster(MMM, xmn, xmx, ymn,  ymx, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#   # plot(r)
#   name_time <- time + TS[i]
#   name_time <- str_pad(name_time, 4, pad="0")
#   name_time <- str_sub(name_time, start = 3, end = -1)
#   name_time <- paste0(name_time, "_00")
#   
#   # crop over the UAE area
#   r <- crop(r, extent(shp_WRF))
#   r <- mask(r, shp_WRF)
#   #   plot(r)
#   
#   # kirging
#   ### Extract points from raster tiff ############################################
#   
#   r_pts <- rasterToPoints(r)
#   colnames(r_pts) <- c("lon", "lat", "O3")
#   r_pts <- as.data.frame(r_pts) 
#   
#   # converta data into a spatialpolygon dataframe------------------------------------
#   
#   r_pts$x <- r_pts$lon
#   r_pts$y <- r_pts$lat
#   
#   coordinates(r_pts) = ~x + y  ## Set spatial coordinates to create a Spatial object:
#   
#   
#   # make a variogram----------------------------------------------------------------
#   
#   vargram_ECMWF <- variogram(O3 ~ 1, r_pts) # calculates sample variogram values 
#   
#   # fite the variogram
#   vargram_ECMWF_fit  <- fit.variogram(vargram_ECMWF, fit.ranges = FALSE, fit.sills = FALSE,
#                                       vgm(1, c("Sph", "Exp", "Mat")), fit.kappa = TRUE)
#   
#   # plot(vargram_ECMWF, vargram_ECMWF_fit) # plot the sample values, along with the fit model
#   
#   e <- extent(shp_WRF)
#   
#   # make a regular empty grid
#   x.range <- as.numeric(c(e@xmin, e@xmax))  # min/max longitude of the interpolation area
#   y.range <- as.numeric(c(e@ymin, e@ymax))  # min/max latitude of the interpolation area
#   
#   ## grid at 10km resolution
#   grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1),
#                      y = seq(from = y.range[1], to = y.range[2], by = 0.1))  # expand points to grid
#   coordinates(grd) <- ~x + y
#   gridded(grd) <- TRUE
#   
#   f.1 <- as.formula(Precip_in ~ X + Y)
#   # perform kriging
#   dat.krg <- gstat::krige(O3 ~ 1, r_pts, grd, vargram_ECMWF_fit, nmax = 50)
#   r <- raster(dat.krg)
#   projection(r) <- CRS("+proj=longlat +datum=WGS84")
#   # r <- crop(r, extent(shp_UAE))
#   # r <- mask(r, shp_UAE)
#   
#   # plot(r)
#   # plot(shp_UAE, add = T, lwd =1)
#   
#   
#   
#   #########-------------------------------------------------------------------------
#   # load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)
#   
#   PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   #   PM25_2015_GRW <- raster("D:/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")
#   
#   ### crop raster over the UAE shp file  ###############################
#   PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
#   PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)
#   
#   ###################
#   
#   ## make resolution of MODIS-data as the one of GWR-------------------------------
#   ECMWF_tif_1km = projectRaster(r, PM25_2015_GRW)
#   ECMWF_tif_1km <- crop(ECMWF_tif_1km, extent(shp_UAE))
#   ECMWF_tif_1km <- mask(ECMWF_tif_1km, shp_UAE)
#   
#   #  plot(ECMWF_tif_1km)
#   
#   #  writeRaster(ECMWF_tif_1km, paste0("D:/ECMWF_data_download/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   writeRaster(ECMWF_tif_1km, paste0(wd,"/",folders[1],"/O3/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   
# }
# 
# 
# 
# #####################################################################################
# #####################################################################################
# ## get NO2 ###########################################################################
# #####################################################################################
# #####################################################################################
# 
# # setwd("D:/ECMWF_data_download")
# setwd("/disk3/fkaragulian/ECMWF_forecasts/forecasts_downloads/")
# 
# 
# # dir <- "D:/MODIS_AOD/UAE_boundary"
# dir <- "/disk3/fkaragulian/MODIS_AOD/UAE_boundary"
# ### shapefile for UAE
# shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")
# 
# # dir <- "D:/Dust_Event_UAE_2015/WRFChem_domain"
# 
# # larger WRF domain
# shp_WRF <- readOGR(dsn = dir, layer = "domain_d02_4km_WRFChem_small")
# # plot(shp_WRF)
# # plot(shp_UAE, add=TRUE, lwd=1)
# 
# 
# 
# folders <- list.files()
# # get the most recent directory
# folders <- sort(folders, decreasing = T)
# 
# wd <- getwd()
# # go to the PM25 folder
# setwd(paste0(wd,"/",folders[1],"/NO2"))
# 
# 
# patt<- ".nc"
# filenames <- list.files(pattern = patt)
# # filenames <- filenames[1]
# 
# 
# # gerate a time sequence of 12 hours
# TS <- 1:12
# 
# # i <- 1
# 
# for(i in 1:length(filenames)) {
#   ECMWF_file <- open.nc(filenames[i])
#   ECMWF_file <- read.nc(ECMWF_file)
#   name_vari <- names(ECMWF_file)
#   
#   
#   
#   datetime <- ECMWF_file[3] 
#   datetime <- as.numeric(datetime)
#   datetime <- convertDateNcdf2R(datetime,units="hours",origin = as.POSIXct("1900-01-01 00:00:01", tz = "UTC"), 
#                                 time.format = c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S", 
#                                                 "%Y-%m-%d %H:%M", "%Y-%m-%d %Z %H:%M", "%Y-%m-%d %Z %H:%M:%S"))
#   datetime <- as.character(datetime)
#   time <- as.numeric(str_sub(filenames[i], start = 23, end = -30))
#   # time <- 0
#   # time <- 12
#   if (time == 0) {
#     time <- 11
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   } else {
#     time <- 0
#     
#     # year <- str_sub(datetime, start = 0, end = -7)
#     # month <- str_sub(datetime, start = 6, end = -4)
#     # day <- str_sub(datetime, start = 9, end = -1) 
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   }
#   
#   
#   NO2 <- (ECMWF_file[4])    #  only one variable (tcno2) == var = 4
#   names(NO2)<- "xxyyzz"
#   NO2 <- (NO2$xxyyzz)
#   LON <- (ECMWF_file$longitude) -180
#   # LON_West <- LON[1:451] -180
#   # LON_East <- LON[452:900] -180
#   # LON <- c(LON_West, LON_East)
#   LAT <- ECMWF_file$latitude
#   
#   
#   xmn= min(LON)
#   xmx=max(LON)
#   ymn=min(LAT)
#   ymx=max(LAT)
#   
#   MMM <-  NO2[ , ]*1e+09
#   MMM <- MMM[ c((451:1),(900:452)), ]
#   MMM <- MMM[nrow(MMM):1, ]
#   MMM <-  t(MMM[ , ])  # map is upside down 
#   r <- raster(MMM, xmn, xmx, ymn,  ymx, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#   # plot(r)
#   name_time <- time + TS[i]
#   name_time <- str_pad(name_time, 4, pad="0")
#   name_time <- str_sub(name_time, start = 3, end = -1)
#   name_time <- paste0(name_time, "_00")
#   
#   # crop over the UAE area
#   r <- crop(r, extent(shp_WRF))
#   r <- mask(r, shp_WRF)
#   #   plot(r)
#   
#   # kirging
#   ### Extract points from raster tiff ############################################
#   
#   r_pts <- rasterToPoints(r)
#   colnames(r_pts) <- c("lon", "lat", "NO2")
#   r_pts <- as.data.frame(r_pts) 
#   
#   # converta data into a spatialpolygon dataframe------------------------------------
#   
#   r_pts$x <- r_pts$lon
#   r_pts$y <- r_pts$lat
#   
#   coordinates(r_pts) = ~x + y  ## Set spatial coordinates to create a Spatial object:
#   
#   
#   # make a variogram----------------------------------------------------------------
#   
#   vargram_ECMWF <- variogram(NO2 ~ 1, r_pts) # calculates sample variogram values 
#   
#   # fite the variogram
#   vargram_ECMWF_fit  <- fit.variogram(vargram_ECMWF, fit.ranges = FALSE, fit.sills = FALSE,
#                                       vgm(1, c("Sph", "Exp", "Mat")), fit.kappa = TRUE)
#   
#   # plot(vargram_ECMWF, vargram_ECMWF_fit) # plot the sample values, along with the fit model
#   
#   e <- extent(shp_WRF)
#   
#   # make a regular empty grid
#   x.range <- as.numeric(c(e@xmin, e@xmax))  # min/max longitude of the interpolation area
#   y.range <- as.numeric(c(e@ymin, e@ymax))  # min/max latitude of the interpolation area
#   
#   ## grid at 10km resolution
#   grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1),
#                      y = seq(from = y.range[1], to = y.range[2], by = 0.1))  # expand points to grid
#   coordinates(grd) <- ~x + y
#   gridded(grd) <- TRUE
#   
#   f.1 <- as.formula(Precip_in ~ X + Y)
#   # perform kriging
#   dat.krg <- gstat::krige(NO2 ~ 1, r_pts, grd, vargram_ECMWF_fit, nmax = 50)
#   r <- raster(dat.krg)
#   projection(r) <- CRS("+proj=longlat +datum=WGS84")
#   # r <- crop(r, extent(shp_UAE))
#   # r <- mask(r, shp_UAE)
#   
#   # plot(r)
#   # plot(shp_UAE, add = T, lwd =1)
#   
#   
#   
#   #########-------------------------------------------------------------------------
#   # load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)
#   
#   PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   #   PM25_2015_GRW <- raster("D:/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")
#   
#   ### crop raster over the UAE shp file  ###############################
#   PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
#   PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)
#   
#   ###################
#   
#   ## make resolution of MODIS-data as the one of GWR-------------------------------
#   ECMWF_tif_1km = projectRaster(r, PM25_2015_GRW)
#   ECMWF_tif_1km <- crop(ECMWF_tif_1km, extent(shp_UAE))
#   ECMWF_tif_1km <- mask(ECMWF_tif_1km, shp_UAE)
#   
#   #  plot(ECMWF_tif_1km)
#   
#   #  writeRaster(ECMWF_tif_1km, paste0("D:/ECMWF_data_download/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   writeRaster(ECMWF_tif_1km, paste0(wd,"/",folders[1],"/NO2/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   
# }
# 
# 
# 
# 
# #####################################################################################
# #####################################################################################
# ## get CO ###########################################################################
# #####################################################################################
# #####################################################################################
# 
# # setwd("D:/ECMWF_data_download")
# setwd("/disk3/fkaragulian/ECMWF_forecasts/forecasts_downloads/")
# 
# 
# # dir <- "D:/MODIS_AOD/UAE_boundary"
# dir <- "/disk3/fkaragulian/MODIS_AOD/UAE_boundary"
# ### shapefile for UAE
# shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")
# 
# # dir <- "D:/Dust_Event_UAE_2015/WRFChem_domain"
# 
# # larger WRF domain
# shp_WRF <- readOGR(dsn = dir, layer = "domain_d02_4km_WRFChem_small")
# # plot(shp_WRF)
# # plot(shp_UAE, add=TRUE, lwd=1)
# 
# 
# 
# folders <- list.files()
# # get the most recent directory
# folders <- sort(folders, decreasing = T)
# 
# wd <- getwd()
# # go to the PM25 folder
# setwd(paste0(wd,"/",folders[1],"/CO"))
# 
# 
# patt<- ".nc"
# filenames <- list.files(pattern = patt)
# # filenames <- filenames[1]
# 
# 
# # gerate a time sequence of 12 hours
# TS <- 1:12
# 
# # i <- 1
# 
# for(i in 1:length(filenames)) {
#   ECMWF_file <- open.nc(filenames[i])
#   ECMWF_file <- read.nc(ECMWF_file)
#   name_vari <- names(ECMWF_file)
#   
#   
#   
#   datetime <- ECMWF_file[3] 
#   datetime <- as.numeric(datetime)
#   datetime <- convertDateNcdf2R(datetime,units="hours",origin = as.POSIXct("1900-01-01 00:00:01", tz = "UTC"), 
#                                 time.format = c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S", 
#                                                 "%Y-%m-%d %H:%M", "%Y-%m-%d %Z %H:%M", "%Y-%m-%d %Z %H:%M:%S"))
#   datetime <- as.character(datetime)
#   time <- as.numeric(str_sub(filenames[i], start = 23, end = -30))
#   # time <- 0
#   # time <- 12
#   if (time == 0) {
#     time <- 11
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   } else {
#     time <- 0
#     
#     # year <- str_sub(datetime, start = 0, end = -7)
#     # month <- str_sub(datetime, start = 6, end = -4)
#     # day <- str_sub(datetime, start = 9, end = -1) 
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   }
#   
#   
#   CO <- (ECMWF_file[4])    #  only one variable (tcco) == var = 4
#   names(CO)<- "xxyyzz"
#   CO <- (CO$xxyyzz)
#   LON <- (ECMWF_file$longitude) -180
#   # LON_West <- LON[1:451] -180
#   # LON_East <- LON[452:900] -180
#   # LON <- c(LON_West, LON_East)
#   LAT <- ECMWF_file$latitude
#   
#   
#   xmn= min(LON)
#   xmx=max(LON)
#   ymn=min(LAT)
#   ymx=max(LAT)
#   
#   MMM <-  CO[ , ]*1e+09
#   MMM <- MMM[ c((451:1),(900:452)), ]
#   MMM <- MMM[nrow(MMM):1, ]
#   MMM <-  t(MMM[ , ])  # map is upside down 
#   r <- raster(MMM, xmn, xmx, ymn,  ymx, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#   # plot(r)
#   name_time <- time + TS[i]
#   name_time <- str_pad(name_time, 4, pad="0")
#   name_time <- str_sub(name_time, start = 3, end = -1)
#   name_time <- paste0(name_time, "_00")
#   
#   # crop over the UAE area
#   r <- crop(r, extent(shp_WRF))
#   r <- mask(r, shp_WRF)
#   #   plot(r)
#   
#   # kirging
#   ### Extract points from raster tiff ############################################
#   
#   r_pts <- rasterToPoints(r)
#   colnames(r_pts) <- c("lon", "lat", "CO")
#   r_pts <- as.data.frame(r_pts) 
#   
#   # converta data into a spatialpolygon dataframe------------------------------------
#   
#   r_pts$x <- r_pts$lon
#   r_pts$y <- r_pts$lat
#   
#   coordinates(r_pts) = ~x + y  ## Set spatial coordinates to create a Spatial object:
#   
#   
#   # make a variogram----------------------------------------------------------------
#   
#   vargram_ECMWF <- variogram(CO ~ 1, r_pts) # calculates sample variogram values 
#   
#   # fite the variogram
#   vargram_ECMWF_fit  <- fit.variogram(vargram_ECMWF, fit.ranges = FALSE, fit.sills = FALSE,
#                                       vgm(1, c("Sph", "Exp", "Mat")), fit.kappa = TRUE)
#   
#   # plot(vargram_ECMWF, vargram_ECMWF_fit) # plot the sample values, along with the fit model
#   
#   e <- extent(shp_WRF)
#   
#   # make a regular empty grid
#   x.range <- as.numeric(c(e@xmin, e@xmax))  # min/max longitude of the interpolation area
#   y.range <- as.numeric(c(e@ymin, e@ymax))  # min/max latitude of the interpolation area
#   
#   ## grid at 10km resolution
#   grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1),
#                      y = seq(from = y.range[1], to = y.range[2], by = 0.1))  # expand points to grid
#   coordinates(grd) <- ~x + y
#   gridded(grd) <- TRUE
#   
#   f.1 <- as.formula(Precip_in ~ X + Y)
#   # perform kriging
#   dat.krg <- gstat::krige(CO ~ 1, r_pts, grd, vargram_ECMWF_fit, nmax = 50)
#   r <- raster(dat.krg)
#   projection(r) <- CRS("+proj=longlat +datum=WGS84")
#   # r <- crop(r, extent(shp_UAE))
#   # r <- mask(r, shp_UAE)
#   
#   # plot(r)
#   # plot(shp_UAE, add = T, lwd =1)
#   
#   
#   
#   #########-------------------------------------------------------------------------
#   # load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)
#   
#   PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   #   PM25_2015_GRW <- raster("D:/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")
#   
#   ### crop raster over the UAE shp file  ###############################
#   PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
#   PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)
#   
#   ###################
#   
#   ## make resolution of MODIS-data as the one of GWR-------------------------------
#   ECMWF_tif_1km = projectRaster(r, PM25_2015_GRW)
#   ECMWF_tif_1km <- crop(ECMWF_tif_1km, extent(shp_UAE))
#   ECMWF_tif_1km <- mask(ECMWF_tif_1km, shp_UAE)
#   
#   #  plot(ECMWF_tif_1km)
#   
#   #  writeRaster(ECMWF_tif_1km, paste0("D:/ECMWF_data_download/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   writeRaster(ECMWF_tif_1km, paste0(wd,"/",folders[1],"/CO/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   
# }
# 
# 
# 
# #####################################################################################
# #####################################################################################
# ## get BLACK CARBON AOD #############################################################
# #####################################################################################
# #####################################################################################
# 
# # setwd("D:/ECMWF_data_download")
# setwd("/disk3/fkaragulian/ECMWF_forecasts/forecasts_downloads/")
# 
# 
# # dir <- "D:/MODIS_AOD/UAE_boundary"
# dir <- "/disk3/fkaragulian/MODIS_AOD/UAE_boundary"
# ### shapefile for UAE
# shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")
# 
# # dir <- "D:/Dust_Event_UAE_2015/WRFChem_domain"
# 
# # larger WRF domain
# shp_WRF <- readOGR(dsn = dir, layer = "domain_d02_4km_WRFChem_small")
# # plot(shp_WRF)
# # plot(shp_UAE, add=TRUE, lwd=1)
# 
# 
# 
# folders <- list.files()
# # get the most recent directory
# folders <- sort(folders, decreasing = T)
# 
# wd <- getwd()
# # go to the PM25 folder
# setwd(paste0(wd,"/",folders[1],"/BC_AOD"))
# 
# 
# patt<- ".nc"
# filenames <- list.files(pattern = patt)
# # filenames <- filenames[1]
# 
# 
# # gerate a time sequence of 12 hours
# TS <- 1:12
# 
# # i <- 1
# 
# for(i in 1:length(filenames)) {
#   ECMWF_file <- open.nc(filenames[i])
#   ECMWF_file <- read.nc(ECMWF_file)
#   name_vari <- names(ECMWF_file)
#   
#   
#   
#   datetime <- ECMWF_file[3] 
#   datetime <- as.numeric(datetime)
#   datetime <- convertDateNcdf2R(datetime,units="hours",origin = as.POSIXct("1900-01-01 00:00:01", tz = "UTC"), 
#                                 time.format = c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S", 
#                                                 "%Y-%m-%d %H:%M", "%Y-%m-%d %Z %H:%M", "%Y-%m-%d %Z %H:%M:%S"))
#   datetime <- as.character(datetime)
#   time <- as.numeric(str_sub(filenames[i], start = 23, end = -30))
#   # time <- 0
#   # time <- 12
#   if (time == 0) {
#     time <- 11
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   } else {
#     time <- 0
#     
#     # year <- str_sub(datetime, start = 0, end = -7)
#     # month <- str_sub(datetime, start = 6, end = -4)
#     # day <- str_sub(datetime, start = 9, end = -1) 
#     
#     year <- str_sub(datetime, start = 0, end = -16)
#     month <- str_sub(datetime, start = 6, end = -13)
#     day <- str_sub(datetime, start = 9, end = -10) 
#     
#     date <- paste0(year,"_",month,"_", day)
#     
#   }
#   
#   
#   BC_AOD <- (ECMWF_file[4])    #  only one variable (tcco) == var = 4
#   names(BC_AOD)<- "xxyyzz"
#   BC_AOD <- (BC_AOD$xxyyzz)
#   LON <- (ECMWF_file$longitude) -180
#   # LON_West <- LON[1:451] -180
#   # LON_East <- LON[452:900] -180
#   # LON <- c(LON_West, LON_East)
#   LAT <- ECMWF_file$latitude
#   
#   
#   xmn= min(LON)
#   xmx=max(LON)
#   ymn=min(LAT)
#   ymx=max(LAT)
#   
#   MMM <-  BC_AOD[ , ]
#   MMM <- MMM[ c((451:1),(900:452)), ]
#   MMM <- MMM[nrow(MMM):1, ]
#   MMM <-  t(MMM[ , ])  # map is upside down 
#   r <- raster(MMM, xmn, xmx, ymn,  ymx, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#   # plot(r)
#   name_time <- time + TS[i]
#   name_time <- str_pad(name_time, 4, pad="0")
#   name_time <- str_sub(name_time, start = 3, end = -1)
#   name_time <- paste0(name_time, "_00")
#   
#   # crop over the UAE area
#   r <- crop(r, extent(shp_WRF))
#   r <- mask(r, shp_WRF)
#   #   plot(r)
#   
#   # kirging
#   ### Extract points from raster tiff ############################################
#   
#   r_pts <- rasterToPoints(r)
#   colnames(r_pts) <- c("lon", "lat", "BC_AOD")
#   r_pts <- as.data.frame(r_pts) 
#   
#   # converta data into a spatialpolygon dataframe------------------------------------
#   
#   r_pts$x <- r_pts$lon
#   r_pts$y <- r_pts$lat
#   
#   coordinates(r_pts) = ~x + y  ## Set spatial coordinates to create a Spatial object:
#   
#   
#   # make a variogram----------------------------------------------------------------
#   
#   vargram_ECMWF <- variogram(BC_AOD ~ 1, r_pts) # calculates sample variogram values 
#   
#   # fite the variogram
#   vargram_ECMWF_fit  <- fit.variogram(vargram_ECMWF, fit.ranges = FALSE, fit.sills = FALSE,
#                                       vgm(1, c("Sph", "Exp", "Mat")), fit.kappa = TRUE)
#   
#   # plot(vargram_ECMWF, vargram_ECMWF_fit) # plot the sample values, along with the fit model
#   
#   e <- extent(shp_WRF)
#   
#   # make a regular empty grid
#   x.range <- as.numeric(c(e@xmin, e@xmax))  # min/max longitude of the interpolation area
#   y.range <- as.numeric(c(e@ymin, e@ymax))  # min/max latitude of the interpolation area
#   
#   ## grid at 10km resolution
#   grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1),
#                      y = seq(from = y.range[1], to = y.range[2], by = 0.1))  # expand points to grid
#   coordinates(grd) <- ~x + y
#   gridded(grd) <- TRUE
#   
#   f.1 <- as.formula(Precip_in ~ X + Y)
#   # perform kriging
#   dat.krg <- gstat::krige(BC_AOD ~ 1, r_pts, grd, vargram_ECMWF_fit, nmax = 50)
#   r <- raster(dat.krg)
#   projection(r) <- CRS("+proj=longlat +datum=WGS84")
#   # r <- crop(r, extent(shp_UAE))
#   # r <- mask(r, shp_UAE)
#   
#   # plot(r)
#   # plot(shp_UAE, add = T, lwd =1)
#   
#   
#   
#   #########-------------------------------------------------------------------------
#   # load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)
#   
#   PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   #   PM25_2015_GRW <- raster("D:/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
#   projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")
#   
#   ### crop raster over the UAE shp file  ###############################
#   PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
#   PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)
#   
#   ###################
#   
#   ## make resolution of MODIS-data as the one of GWR-------------------------------
#   ECMWF_tif_1km = projectRaster(r, PM25_2015_GRW)
#   ECMWF_tif_1km <- crop(ECMWF_tif_1km, extent(shp_UAE))
#   ECMWF_tif_1km <- mask(ECMWF_tif_1km, shp_UAE)
#   
#   #  plot(ECMWF_tif_1km)
#   
#   #  writeRaster(ECMWF_tif_1km, paste0("D:/ECMWF_data_download/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   writeRaster(ECMWF_tif_1km, paste0(wd,"/",folders[1],"/BC_AOD/", date, "_", name_time, ".tif") , options= "INTERLEAVE=BAND", overwrite=T)
#   
# }


