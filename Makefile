ecoregion_summaries.csv: R/aggregate-climate-data.R R/get-ecoregion-summaries.R
	Rscript --vanilla R/aggregate-climate-data.R
	Rscript --vanilla R/get-ecoregion-summaries.R
	aws s3 cp ecoregion_summaries.csv s3://earthlab-gridmet/ecoregion_summaries.csv --acl public-read
