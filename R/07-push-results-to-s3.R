
# Push important stuff to S3 ----------------------------------------------
model_fits <- list.files(pattern = '*.fit.*\\.rds')

for (i in seq_along(model_fits)) {
  cmd <- paste0('aws s3 cp ',
               model_fits[i],
               ' s3://earthlab-mjoseph/demo_evt/',
               model_fits[i])
  system(cmd)
}

# push figures to S3
#system('aws s3 cp fig s3://earthlab-mjoseph/demo_evt/fig --recursive')

# pull fits down
# system('aws s3 cp s3://earthlab-mjoseph/demo_evt/ba_gamma_fit.rds . ')
# system('aws s3 cp s3://earthlab-mjoseph/demo_evt/ba_lognormal_fit.rds . ')
# system('aws s3 cp s3://earthlab-mjoseph/demo_evt/ba_pareto_fit.rds . ')
# system('aws s3 cp s3://earthlab-mjoseph/demo_evt/ba_tpareto_fit.rds . ')
# system('aws s3 cp s3://earthlab-mjoseph/demo_evt/ba_weibull_fit.rds . ')

# system('aws s3 cp s3://earthlab-mjoseph/demo_evt/pois_fit.rds . ')
# system('aws s3 cp s3://earthlab-mjoseph/demo_evt/nb_fit.rds . ')

