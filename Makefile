
data: data/processed/ecoregion_summaries.csv data/processed/housing_density.csv data/processed/ecoregion_tri.csv

data/processed/ecoregion_summaries.csv: R/aggregate-climate-data.R R/get-ecoregion-summaries.R R/fetch-fire-data.R
	Rscript --vanilla R/fetch-fire-data.R
	Rscript --vanilla R/aggregate-climate-data.R
	Rscript --vanilla R/get-ecoregion-summaries.R
	aws s3 cp data/processed/ecoregion_summaries.csv s3://earthlab-gridmet/ecoregion_summaries.csv --acl public-read

data/processed/housing_density.csv: R/summarize-housing-density.R R/helpers.R
	# download the data from the internet:
	wget -O data/raw/us_pbg00.zip http://silvis.forest.wisc.edu/sites/default/files/maps/pbg00_old/gis/us_pbg00.zip
	unzip data/raw/us_pbg00.zip -d data/raw/

	# use GDAL to convert to an ESRI Shapefile
	ogr2ogr -f 'ESRI Shapefile' data/processed/us_pbg data/raw/us_pbg00_2007.gdb

	# rasterize for each year of interest 1980 - 2019 (4km grid):
	SHP="data/processed/us_pbg/us_pbg00_Oct_8_2007.shp"
	gdal_rasterize -a HDEN80 -of GTiff -tr 4000 4000 $SHP data/processed/den80.tif
	gdal_rasterize -a HDEN90 -of GTiff -tr 4000 4000 $SHP data/processed/den90.tif
	gdal_rasterize -a HDEN00 -of GTiff -tr 4000 4000 $SHP data/processed/den00.tif
	gdal_rasterize -a HDEN10 -of GTiff -tr 4000 4000 $SHP data/processed/den10.tif
	gdal_rasterize -a HDEN20 -of GTiff -tr 4000 4000 $SHP data/processed/den20.tif

  Rscript --vanilla R/summarize-housing-density.R
	aws s3 cp data/processed/housing_density.csv s3://earthlab-mjoseph/demo_evt/housing_density.csv --acl public-read

data/processed/ecoregion_tri.csv: R/compute-terrain-metrics.R
	Rscript --vanilla R/compute-terrain-metrics.R
