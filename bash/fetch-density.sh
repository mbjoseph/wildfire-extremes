#!/bin/bash


# Download the data from the internet:

wget -O data/raw/us_pbg00.zip http://silvis.forest.wisc.edu/sites/default/files/maps/pbg00_old/gis/us_pbg00.zip
unzip data/raw/us_pbg00.zip -d data/raw/


# use GDAL to convert to an ESRI Shapefile
ogr2ogr -f 'ESRI Shapefile' data/processed/us_pbg data/raw/us_pbg00_2007.gdb


# rasterize for each year of interest 1980 - 2030 (4km grid):
gdal_rasterize -a HDEN80 -of GTiff -tr 4000 4000 data/processed/us_pbg/us_pbg00_Oct_8_2007.shp data/processed/den80.tif
gdal_rasterize -a HDEN90 -of GTiff -tr 4000 4000 data/processed/us_pbg/us_pbg00_Oct_8_2007.shp data/processed/den90.tif
gdal_rasterize -a HDEN00 -of GTiff -tr 4000 4000 data/processed/us_pbg/us_pbg00_Oct_8_2007.shp data/processed/den00.tif
gdal_rasterize -a HDEN10 -of GTiff -tr 4000 4000 data/processed/us_pbg/us_pbg00_Oct_8_2007.shp data/processed/den10.tif
gdal_rasterize -a HDEN20 -of GTiff -tr 4000 4000 data/processed/us_pbg/us_pbg00_Oct_8_2007.shp data/processed/den20.tif
gdal_rasterize -a HDEN30 -of GTiff -tr 4000 4000 data/processed/us_pbg/us_pbg00_Oct_8_2007.shp data/processed/den30.tif
