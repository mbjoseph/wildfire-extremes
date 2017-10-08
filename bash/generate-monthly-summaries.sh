#!/bin/bash

# download the data from remote server
source metdata_wget.sh

# execute R script to produce monthly GeoTIFFs
Rscript --vanilla R/aggregate-climate-data.R

# push monthly GeoTIFFS to AWS S3
aws s3 cp ../data/processed/rmin s3://earthlab-gridmet/rmin --recursive
aws s3 cp ../data/processed/tmmx s3://earthlab-gridmet/tmmx --recursive
aws s3 cp ../data/processed/vs s3://earthlab-gridmet/vs --recursive
aws s3 cp ../data/processed/pr s3://earthlab-gridmet/pr --recursive

