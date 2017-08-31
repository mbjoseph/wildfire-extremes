# Spatiotemporal modeling of wildfire extremes

## Getting the data

1. `bash/generate-monthly-summaries.sh` downloads climate data from MET, aggregates to a monthly timestep, and pushes the resulting GeoTIFF files to Amazon S3. This step requires a fair bit of memory (64 GB).

2. `bash/generate-ecoregion-summaries.sh` downloads monthly summaries from S3, computes US EPA level 3 ecoregion summaries, and pushes the resulting tidy csv file back to S3. This step requires 40 CPU cores. 

3. `R/fetch-fire-data.R` downloads EPA level 3 ecoregion and Monitoring Trends in Burn Severity (MTBS) fire shapefiles. 

## Cleaning data

4. `R/01-clean_data.R` will load the fire data and merge it with the ecoregion and climate data.

