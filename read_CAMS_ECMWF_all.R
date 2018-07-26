
library(RNetCDF)
library(fields)
library(maptools)
library(rgdal)
library(readr)
library(ggplot2)
library(gstat)
library(sp)
library(maptools)
library(raster)   
library(rgeos)
library(plyr)
library(dplyr)
library(leaflet)
library(htmltools)
library(ncdf4)
library(htmlwidgets)
library(spatialEco)
library(pracma)
library(threadr)
library(webshot)
library(readr)

setwd ("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2015") 

dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/website_MODIS/UAE_boundary"
### shapefile for UAE
shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")

# ----- Transform to EPSG 4326 - WGS84 (required)
shp_UAE <- spTransform(shp_UAE, CRS("+init=epsg:4326"))
# names(shp)

shp_UAE@data$name <- 1:nrow(shp_UAE)
# plot(shp_UAE)

filenames <- list.files(pattern = "\\.nc") # list all forecast for next 5 days (ECMWF)
# file <- filenames[1]
# i <- 1

# FUNCTION to read data for each month & day---------------------------------------

ECMWF_DAY <- NULL


ECMWF_DAYS <- function (file) {      ## this is the filenames 
  date <- str_sub(file, start = 1, end = -4)
#  date <- str_sub(filenames[1], start = 1, end = -4)
  ECMWF_file <- file

  # open raster file and crop it  
  ECMWF_file <- raster(file)
  ECMWF_file <- crop(ECMWF_file, extent(shp_UAE))
  ECMWF_file <- mask(ECMWF_file, shp_UAE)
  plot(ECMWF_file)
  res(ECMWF_file)
  
  
  
#  ECMWF_file <- open.nc(filenames[2])
ECMWF_file <- open.nc(file)
ECMWF_file <- read.nc(ECMWF_file)

# exctract the data for PM2.5-------------------------------------------
LAT_ECMWF <- as.vector(ECMWF_file$latitude)  ### get Lat
LON_ECMWF <- as.vector(ECMWF_file$longitude)  ### get Lon
LAT <- repmat(LAT_ECMWF,1,26)
LAT <- t(LAT)
LAT <- as.data.frame(LAT)
LON = repmat(LON_ECMWF,21,1)
LON = c(LON)
LON <- as.data.frame(LON)
colnames(LON) <- "LON"
colnames(LAT) <- "LAT"

# DATA = NULL

## loop thorugh all days for PM2.5 --------------
for (i in 1: ncol(ECMWF_file$pm2p5[1,,])) {
  
AOD <- ECMWF_file$aod550[,,i]  ### AOD 
PM25 <- (ECMWF_file$pm2p5[,,i])*1000000000 ### PM2.5 (from kg/m3 to ug/m3)
PM10 <- (ECMWF_file$pm10[,,i])*1000000000 ### PM2.5 (from kg/m3 to ug/m3)
DUST_AOD <- (ECMWF_file$duaod550[,,i]) ### DUST (AOD)
SALT_AOD <- (ECMWF_file$ssaod550[,,i]) ### SEA SALT (AOD)
BC_AOD <- (ECMWF_file$bcaod550[,,i]) ### Black Carbon (AOD)
OM_AOD <- (ECMWF_file$bcaod550[,,i]) ### Organic Matter (AOD)
SULPHATE_AOD <- (ECMWF_file$suaod550[,,i]) ### Organic Matter (AOD)

AOD <- t(AOD)
AOD <- c(AOD)
AOD <- as.data.frame(AOD)

PM10 <- t(PM10)
PM10 <- c(PM10)
PM10 <- as.data.frame(PM10)

PM25 <- t(PM25)
PM25 <- c(PM25)
PM25 <- as.data.frame(PM25)

DUST_AOD <- t(DUST_AOD)
DUST_AOD <- c(DUST_AOD)
DUST_AOD <- as.data.frame(DUST_AOD)

SALT_AOD <- t(SALT_AOD)
SALT_AOD <- c(SALT_AOD)
SALT_AOD <- as.data.frame(SALT_AOD)

BC_AOD <- t(BC_AOD)
BC_AOD <- c(BC_AOD)
BC_AOD <- as.data.frame(BC_AOD)

OM_AOD <- t(OM_AOD)
OM_AOD <- c(OM_AOD)
OM_AOD <- as.data.frame(OM_AOD)

SULPHATE_AOD <- t(SULPHATE_AOD)
SULPHATE_AOD <- c(SULPHATE_AOD)
SULPHATE_AOD <- as.data.frame(SULPHATE_AOD)


# bind lon, lat, data
ECMWF_day <- cbind(LON, LAT, AOD, PM25, PM10, DUST_AOD, SALT_AOD, BC_AOD, OM_AOD, SULPHATE_AOD)

# make an interpolation not a regular grid---------
ECMWF_day$x <- ECMWF_day$LON
ECMWF_day$y <- ECMWF_day$LAT

coordinates(ECMWF_day) = ~ x + y  ## Set spatial coordinates to create a Spatial object:


# make a regular empty grid
x.range <- as.numeric(c(20, 70))  # min/max longitude of the interpolation area
y.range <- as.numeric(c(-2, 50))  # min/max latitude of the interpolation area


## grid at 10km resolution
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.4),
                   y = seq(from = y.range[1], to = y.range[2], by = 0.4))  # expand points to grid
coordinates(grd) <- ~x + y
gridded(grd) <- TRUE

plot(grd, cex = 1.5, col = "grey")
points(ECMWF_day, pch = 1, col = "red", cex = 1)


idw_AOD <- idw(formula = AOD ~ 1, locations = ECMWF_day,
                newdata = grd)  # apply idw model for the data (interpolation)

idw_PM25 <- idw(formula = PM25 ~ 1, locations = ECMWF_day,
           newdata = grd)  

idw_PM10 <- idw(formula = PM10 ~ 1, locations = ECMWF_day,
                newdata = grd)  

idw_DUST_AOD <- idw(formula = DUST_AOD ~ 1, locations = ECMWF_day,
               newdata = grd)  

idw_SALT_AOD <- idw(formula = SALT_AOD ~ 1, locations = ECMWF_day,
                    newdata = grd)  

idw_BC_AOD <- idw(formula = BC_AOD ~ 1, locations = ECMWF_day,
                    newdata = grd) 

idw_OM_AOD <- idw(formula = OM_AOD ~ 1, locations = ECMWF_day,
                  newdata = grd) 

idw_SULPHATE_AOD <- idw(formula = SULPHATE_AOD ~ 1, locations = ECMWF_day,
                  newdata = grd) 

idw.output_AOD = as.data.frame(idw_AOD)  
idw.output_PM25 = as.data.frame(idw_PM25)
idw.output_PM10 = as.data.frame(idw_PM10)
idw.output_DUST_AOD = as.data.frame(idw_DUST_AOD)
idw.output_SALT_AOD = as.data.frame(idw_SALT_AOD)
idw.output_BC_AOD = as.data.frame(idw_BC_AOD)
idw.output_OM_AOD = as.data.frame(idw_OM_AOD)
idw.output_SULPHATE_AOD = as.data.frame(idw_SULPHATE_AOD)

names(idw.output_AOD)[1:3] <- c("Lon", "Lat", "AOD")  # give names to the modelled variables
names(idw.output_PM25)[1:3] <- c("Lon", "Lat", "PM25")  
names(idw.output_PM10)[1:3] <- c("Lon", "Lat", "PM10")  
names(idw.output_DUST_AOD)[1:3] <- c("Lon", "Lat", "DUST_AOD") 
names(idw.output_SALT_AOD)[1:3] <- c("Lon", "Lat", "SALT_AOD")  
names(idw.output_BC_AOD)[1:3] <- c("Lon", "Lat", "BC_AOD")
names(idw.output_OM_AOD)[1:3] <- c("Lon", "Lat", "OM_AOD")
names(idw.output_SULPHATE_AOD)[1:3] <- c("Lon", "Lat", "SO4_AOD")  



##### make raster AOD----------------------------------------

coordinates(idw.output_AOD) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(idw.output_AOD) <- TRUE
r <- raster(idw.output_AOD)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
r <- crop(r, extent(shp_UAE))
r <- mask(r, shp_UAE)
plot(r)


# extract points

ECMWF_day_pts_AOD <- rasterToPoints(r)
head(ECMWF_day_pts_AOD)
colnames(ECMWF_day_pts_AOD) <- c("Lon", "Lat", "AOD_ECMWF")
ECMWF_day_pts_AOD <- as.data.frame (ECMWF_day_pts_AOD)
# add column for the date
ECMWF_day_pts_AOD$Date <- paste0(date,"-", sprintf("%02d", i))


# # save csv file for each day-----------------------
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/AOD_ECMWF")
# dir_AOD <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/AOD_ECMWF"
# write_csv(ECMWF_day_pts_AOD, paste0(dir_AOD, "/", date,"-", sprintf("%02d", i), ".csv"))


#### to make raster AOD at 1km------------------------------------

# coordinates(idw.output_AOD) <- ~ Lon + Lat
# # coerce to SpatialPixelsDataFrame
# gridded(idw.output_AOD) <- TRUE
# r <- raster(idw.output_AOD)
# projection(r) <- CRS("+proj=longlat +datum=WGS84")
# 
# # load a MODIS raster at 1km resolution  
# MODIS_1km <- raster("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed/2016_AOD_tiff_1km/2016-01-01.tif")
# # res(MODIS_1km)
# # plot(MODIS_1km)
# 
# # regridd idw file at 1km using the above raster
# r_1km = projectRaster(r, MODIS_1km) 
# # res(r_1km)
# # plot(r_1km)
# 
# r_1km <- crop(r_1km, extent(shp_UAE))
# r_1km <- mask(r_1km, shp_UAE)
# plot(r_1km)
# 
# 
# # extract points
# 
# ECMWF_day_pts <- rasterToPoints(r_1km)
# head(ECMWF_day_pts)
# colnames(ECMWF_day_pts) <- c("Lon", "Lat", "AOD_ECMWF")
# ECMWF_day_pts <- as.data.frame (ECMWF_day_pts)
# # add column for the date
# ECMWF_day_pts$Date <- paste0(date,"-", sprintf("%02d", i))
# 
# # save csv file for each day----------------------- 
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/AOD_ECMWF")
# dir_AOD <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/AOD_ECMWF"
# write_csv(ECMWF_day_pts, paste0(dir_AOD, "/", date,"-", sprintf("%02d", i), ".csv"))


##### make raster PM25 ----------------------------------------

coordinates(idw.output_PM25) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(idw.output_PM25) <- TRUE
r <- raster(idw.output_PM25)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
r <- crop(r, extent(shp_UAE))
r <- mask(r, shp_UAE)
plot(r)


# extract points

ECMWF_day_pts_PM25 <- rasterToPoints(r)
head(ECMWF_day_pts_PM25)
colnames(ECMWF_day_pts_PM25) <- c("Lon", "Lat", "PM25_ECMWF")
ECMWF_day_pts_PM25 <- as.data.frame (ECMWF_day_pts_PM25)
# add column for the date
ECMWF_day_pts_PM25$Date <- paste0(date,"-", sprintf("%02d", i))

# bind with previous data
ECMWF_DAY <- cbind(ECMWF_day_pts_AOD, ECMWF_day_pts_PM25[3]) # only data


# # save csv file for each day-----------------------
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/PM25_ECMWF")
# dir_PM25 <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/PM25_ECMWF"
# write_csv(ECMWF_day_pts_PM25, paste0(dir_PM25, "/", date,"-", sprintf("%02d", i), ".csv"))


##### make raster PM10 ----------------------------------------

coordinates(idw.output_PM10) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(idw.output_PM10) <- TRUE
r <- raster(idw.output_PM10)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
r <- crop(r, extent(shp_UAE))
r <- mask(r, shp_UAE)
plot(r)


# extract points

ECMWF_day_pts_PM10 <- rasterToPoints(r)
head(ECMWF_day_pts_PM10)
colnames(ECMWF_day_pts_PM10) <- c("Lon", "Lat", "PM10_ECMWF")
ECMWF_day_pts_PM10 <- as.data.frame (ECMWF_day_pts_PM10)
# add column for the date
ECMWF_day_pts_PM10$Date <- paste0(date,"-", sprintf("%02d", i))


# bind with previous data
ECMWF_DAY <- cbind(ECMWF_DAY, ECMWF_day_pts_PM10[3])

# # save csv file for each day-----------------------
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/PM10_ECMWF")
# dir_PM10 <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/PM10_ECMWF"
# write_csv(ECMWF_day_pts_PM10, paste0(dir_PM10, "/", date,"-", sprintf("%02d", i), ".csv"))


##### make raster DUST AOD ----------------------------------------

coordinates(idw.output_DUST_AOD) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(idw.output_DUST_AOD) <- TRUE
r <- raster(idw.output_DUST_AOD)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
r <- crop(r, extent(shp_UAE))
r <- mask(r, shp_UAE)
plot(r)


# extract points

ECMWF_day_pts_DUST <- rasterToPoints(r)
head(ECMWF_day_pts_DUST)
colnames(ECMWF_day_pts_DUST) <- c("Lon", "Lat", "DUST_ECMWF")
ECMWF_day_pts_DUST <- as.data.frame (ECMWF_day_pts_DUST)
# add column for the date
ECMWF_day_pts_DUST$Date <- paste0(date,"-", sprintf("%02d", i))


# bind with previous data
ECMWF_DAY <- cbind(ECMWF_DAY, ECMWF_day_pts_DUST[3])

# # save csv file for each day-----------------------
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/DUST_ECMWF")
# dir_DUST <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/DUST_ECMWF"
# write_csv(ECMWF_day_pts_DUST, paste0(dir_DUST, "/", date,"-", sprintf("%02d", i), ".csv"))


##### make raster SEA SALT AOD ----------------------------------------

coordinates(idw.output_SALT_AOD) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(idw.output_SALT_AOD) <- TRUE
r <- raster(idw.output_SALT_AOD)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
r <- crop(r, extent(shp_UAE))
r <- mask(r, shp_UAE)
plot(r)


# extract points

ECMWF_day_pts_SALT <- rasterToPoints(r)
head(ECMWF_day_pts_SALT)
colnames(ECMWF_day_pts_SALT) <- c("Lon", "Lat", "SALT_ECMWF")
ECMWF_day_pts_SALT <- as.data.frame (ECMWF_day_pts_SALT)
# add column for the date
ECMWF_day_pts_SALT$Date <- paste0(date,"-", sprintf("%02d", i))

# bind with previous data
ECMWF_DAY <- cbind(ECMWF_DAY, ECMWF_day_pts_SALT[3])


# # save csv file for each day-----------------------
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/SALT_ECMWF")
# dir_SALT <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/SALT_ECMWF"
# write_csv(ECMWF_day_pts_SALT, paste0(dir_SALT, "/", date,"-", sprintf("%02d", i), ".csv"))


##### make raster BLACK CARBON AOD ----------------------------------------

coordinates(idw.output_BC_AOD) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(idw.output_BC_AOD) <- TRUE
r <- raster(idw.output_BC_AOD)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
r <- crop(r, extent(shp_UAE))
r <- mask(r, shp_UAE)
plot(r)


# extract points

ECMWF_day_pts_BC <- rasterToPoints(r)
head(ECMWF_day_pts_BC)
colnames(ECMWF_day_pts_BC) <- c("Lon", "Lat", "BC_ECMWF")
ECMWF_day_pts_BC <- as.data.frame (ECMWF_day_pts_BC)
# add column for the date
ECMWF_day_pts_BC$Date <- paste0(date,"-", sprintf("%02d", i))


# bind with previous data
ECMWF_DAY <- cbind(ECMWF_DAY, ECMWF_day_pts_BC[3])

# # save csv file for each day-----------------------
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/BC_ECMWF")
# dir_BC <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/BC_ECMWF"
# write_csv(ECMWF_day_pts_BC, paste0(dir_BC, "/", date,"-", sprintf("%02d", i), ".csv"))


##### make raster ORGANIC MATTER AOD ----------------------------------------

coordinates(idw.output_OM_AOD) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(idw.output_OM_AOD) <- TRUE
r <- raster(idw.output_OM_AOD)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
r <- crop(r, extent(shp_UAE))
r <- mask(r, shp_UAE)
plot(r)


# extract points

ECMWF_day_pts_OM <- rasterToPoints(r)
head(ECMWF_day_pts_OM)
colnames(ECMWF_day_pts_OM) <- c("Lon", "Lat", "OM_ECMWF")
ECMWF_day_pts_OM <- as.data.frame (ECMWF_day_pts_OM)
# add column for the date
ECMWF_day_pts_OM$Date <- paste0(date,"-", sprintf("%02d", i))


# bind with previous data
ECMWF_DAY <- cbind(ECMWF_DAY, ECMWF_day_pts_OM[3])

# # save csv file for each day-----------------------
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/OM_ECMWF")
# dir_OM <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/OM_ECMWF"
# write_csv(ECMWF_day_pts_OM, paste0(dir_OM, "/", date,"-", sprintf("%02d", i), ".csv"))


##### make raster SULPHATE (SO4) AOD ----------------------------------------

coordinates(idw.output_SULPHATE_AOD) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(idw.output_SULPHATE_AOD) <- TRUE
r <- raster(idw.output_SULPHATE_AOD)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
r <- crop(r, extent(shp_UAE))
r <- mask(r, shp_UAE)
plot(r)


# extract points

ECMWF_day_pts_SO4 <- rasterToPoints(r)
head(ECMWF_day_pts_SO4)
colnames(ECMWF_day_pts_SO4) <- c("Lon", "Lat", "SO4_ECMWF")
ECMWF_day_pts_SO4 <- as.data.frame(ECMWF_day_pts_SO4)
# add column for the date
ECMWF_day_pts_SO4$Date <- paste0(date,"-", sprintf("%02d", i))

# bind with previous data
ECMWF_DAY <- cbind(ECMWF_DAY, ECMWF_day_pts_SO4[3])

# # save csv file for each day-----------------------
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/SO4_ECMWF")
# dir_SO4 <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2016/SO4_ECMWF"
# write_csv(ECMWF_day_pts_SO4, paste0(dir_SO4, "/", date,"-", sprintf("%02d", i), ".csv"))

####################################################################################
####################################################################################
## save all data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ECMWF_DAY <- ECMWF_DAY %>%
  dplyr::select(Date,
                Lon,
                Lat,
                AOD_ECMWF,
                PM25_ECMWF,
                PM10_ECMWF,
                DUST_ECMWF,
                SALT_ECMWF,
                BC_ECMWF,
                OM_ECMWF,
                SO4_ECMWF)


dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2015/all_ECMWF")
dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/ECMWF_CAMS/2015/all_ECMWF"
write_csv(ECMWF_DAY, paste0(dir, "/", date, "-", sprintf("%02d", i), ".csv"))

}
}

## bind column data together---------------------
# DATA = cbind(DATA, Data)
## make average of all data (by month)-----------
# mean <- rowMeans(DATA)

BBB <- lapply(filenames, ECMWF_DAYS)




