#!/bin/bash

# pull monthly summaries from S3
aws s3 cp s3://earthlab-gridmet/rmin data/processed/rmin --recursive
aws s3 cp s3://earthlab-gridmet/pr data/processed/pr --recursive
aws s3 cp s3://earthlab-gridmet/tmmx data/processed/tmmx --recursive
aws s3 cp s3://earthlab-gridmet/vs  data/processed/vs --recursive


# summarise at the ecoregion level (in parallel, requires 40 cores)
Rscript --vanilla R/get-ecoregion-summaries.R


# push tidy ecoregion summary file back to s3
aws s3 cp data/processed/ecoregion_summaries.csv s3://earthlab-gridmet/ecoregion_summaries.csv --acl public-read
