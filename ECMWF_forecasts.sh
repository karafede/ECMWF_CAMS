#!/bin/bash


mainPath=/home/fkaragulian/ECMWF_forecasts/ ; dest=/research/cesam/AirQuality/forecasts_ECMWF/
source /home/fkaragulian/export.sh

date=`date +%Y%m%d` ; hr=$1 ; date_1=${date}${hr}
#date=$2 ; hr=$1 ; date_1=${date}${hr}
cd $mainPath

#module load R/3.3.2
#module load netcdf/4.3.2
#module load gdal/2.1.2

# /apps/R/R-3.3.2/bin/Rscript /home/fkaragulian/ECMWF_forecasts/CAMS_NRT_download.R
/apps/R/R-3.3.2/bin/Rscript /home/fkaragulian/ECMWF_forecasts/CAMS_NRT_download_NEW.R


mkdir -p ${dest}/${date_1}/ ; cd ${dest}/${date_1}/ ; 

cp -r ${mainPath}/forecasts_downloads/${date_1}/PM25/*.tif .

scp ${mainPath}/forecasts_downloads/${date_1}/PM25/*.tif fkaragulian@10.102.14.39:/home/fkaragulian/mtg/_SHARED_FOLDERS/AQ_Website/ECMWF_forecasts_2017/PM25/

echo "fine"

exit
