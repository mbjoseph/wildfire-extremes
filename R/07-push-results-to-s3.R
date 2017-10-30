
# Push important stuff to S3 ----------------------------------------------
model_fits <- list.files(pattern = '*.fit.*\\.rds')

for (i in seq_along(model_fits)) {
  cmd <- paste0('aws s3 cp ',
               model_fits[i],
               ' s3://earthlab-mjoseph/demo_evt/',
               model_fits[i])
  system(cmd)
}

system('aws s3 cp fig s3://earthlab-mjoseph/demo_evt/fig --recursive')
system('aws s3 cp loo-ic-table.csv s3://earthlab-mjoseph/demo_evt/loo-ic-table.csv')
