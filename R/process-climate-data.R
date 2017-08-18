system("bash/metdata_wget.sh")
source("R/aggregate-climate-data.R")
system("bash/push-monthly-summaries-to-s3.sh")
#source("R/get-ecoregion-summaries.R")
