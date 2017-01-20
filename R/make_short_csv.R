library(RSQLite)
library(tidyverse)

# load the Short data and save as a csv
# downloaded from https://www.fs.usda.gov/rds/archive/Product/RDS-2013-0009.3/

f <- "data/raw/RDS-2013-0009.3_SQLite/Data/FPA_FOD_20150323.sqlite"
mydb <- dbConnect(SQLite(), f)
dbListTables(mydb)
dbListFields(mydb, "fires")
d <- dbReadTable(mydb, "fires") %>%
  tbl_df()

d %>%
  select(-shape) %>%
  write_csv(path = "data/processed/short-data.csv")
