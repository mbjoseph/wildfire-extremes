S3BUCKET=earthlab-mjoseph

data-dir=data/processed
gdb=data/raw/us_pbg00_2007.gdb

all: data

data: $(data-dir)/ecoregion_summaries.csv $(data-dir)/housing_density.csv

data/processed/ecoregion_summaries.csv: R/aggregate-climate-data.R R/get-ecoregion-summaries.R data/raw/us_eco_l3
	Rscript --vanilla R/aggregate-climate-data.R
	Rscript --vanilla R/get-ecoregion-summaries.R

data/raw/us_pbg00_2007.gdb:
	wget -nc -O data/raw/us_pbg00.zip http://silvis.forest.wisc.edu/sites/default/files/maps/pbg00_old/gis/us_pbg00.zip
	unzip data/raw/us_pbg00.zip -d data/raw/

data/processed/housing_density.csv: R/summarize-housing-density.R data/raw/us_eco_l3 data/raw/us_pbg00_2007.gdb
	gdal_rasterize -a HDEN80 -of GTiff -tr 4000 4000 $(gdb) data/processed/den80.tif
	gdal_rasterize -a HDEN90 -of GTiff -tr 4000 4000 $(gdb) data/processed/den90.tif
	gdal_rasterize -a HDEN00 -of GTiff -tr 4000 4000 $(gdb) data/processed/den00.tif
	gdal_rasterize -a HDEN10 -of GTiff -tr 4000 4000 $(gdb) data/processed/den10.tif
	gdal_rasterize -a HDEN20 -of GTiff -tr 4000 4000 $(gdb) data/processed/den20.tif
	Rscript --vanilla R/summarize-housing-density.R

data/raw/us_eco_l3: R/fetch-fire-data.R
	Rscript --vanilla R/fetch-fire-data.R

s3data: $(processed-data)
	aws s3 cp $(data-dir)/ecoregion_summaries.csv s3://$(S3BUCKET)/ecoregion_summaries.csv --acl public-read
	aws s3 cp $(data-dir)/housing_density.csv s3://$(S3BUCKET)/demo_evt/housing_density.csv --acl public-read
