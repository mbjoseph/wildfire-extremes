data-dir = data/processed
gdb = data/raw/us_pbg00_2007.gdb
figs = fig/ppc-density-funs.png \
	fig/burn-area-effs.png fig/ppc-counts.png \
	fig/all-coefs.png \
	fig/count-partial-effs.png fig/attribution-plot.png \
	fig/count-preds.png \
	fig/max-preds-l3.png \
	fig/maps.png \
	fig/spline-concept.png
	
tables = data/processed/burn-area-loglik.csv data/processed/count-loglik.csv \
	data/processed/rho_beta.csv data/processed/count_test_intervals.csv \
	data/processed/area_df.csv data/processed/predicted_totals.csv \
	data/processed/area_coverage.csv data/processed/mev_intervals.csv \
	data/processed/burn-area-beta.csv

all: main.pdf

main.pdf: $(figs) $(tables) main.Rmd library.bib header.sty
		Rscript -e "rmarkdown::render('main.Rmd')"

data/processed/stan_d.rds data/processed/mtbs.rds data/processed/ecoregions.rds: R/make-stan-d.R \
	$(data-dir)/ecoregion_summaries.csv \
	$(data-dir)/housing_density.csv \
	data/raw/mtbs_fod_pts_data/mtbs_fod_pts_20170501.shp
		Rscript --vanilla R/make-stan-d.R

data/processed/ecoregion_summaries.csv: R/get-ecoregion-summaries.R  \
	data/raw/us_eco_l3/us_eco_l3.shp \
	R/aggregate-climate-data.R \
	data/raw/climate-data.csv
		Rscript --vanilla R/aggregate-climate-data.R
		Rscript --vanilla R/get-ecoregion-summaries.R

data/processed/housing_density.csv: R/summarize-housing-density.R \
	data/raw/us_eco_l3/us_eco_l3.shp \
	data/raw/us_pbg00_2007.gdb
		gdal_rasterize -a HDEN80 -of GTiff -tr 4000 4000 $(gdb) data/processed/den80.tif
		gdal_rasterize -a HDEN90 -of GTiff -tr 4000 4000 $(gdb) data/processed/den90.tif
		gdal_rasterize -a HDEN00 -of GTiff -tr 4000 4000 $(gdb) data/processed/den00.tif
		gdal_rasterize -a HDEN10 -of GTiff -tr 4000 4000 $(gdb) data/processed/den10.tif
		gdal_rasterize -a HDEN20 -of GTiff -tr 4000 4000 $(gdb) data/processed/den20.tif
		Rscript --vanilla R/summarize-housing-density.R

data/raw/us_pbg00_2007.gdb:
		wget -q -nc -O data/raw/us_pbg00.zip http://silvis.forest.wisc.edu/sites/default/files/maps/pbg00_old/gis/us_pbg00.zip
		unzip data/raw/us_pbg00.zip -d data/raw/
		rm data/raw/us_pbg00.zip

data/raw/us_eco_l3/us_eco_l3.shp: 
		wget -q -nc -O data/raw/us_eco_l3.zip ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip
		unzip -o data/raw/us_eco_l3.zip -d data/raw/us_eco_l3
		rm data/raw/us_eco_l3.zip

data/raw/mtbs_fod_pts_data/mtbs_fod_pts_20170501.shp:
		wget -q -nc -O data/raw/mtbs_fod_pts_data.zip https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/fod_pt_shapefile/mtbs_fod_pts_data.zip
		unzip -o data/raw/mtbs_fod_pts_data.zip -d data/raw/mtbs_fod_pts_data
		rm data/raw/mtbs_fod_pts_data.zip

pois_fit.rds: R/fit-count-poisson.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-count-poisson.R

nb_fit.rds: R/fit-count-negbinom.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-count-negbinom.R

zip_fit.rds: R/fit-count-zip.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-count-zip.R

zinb_fit.rds: R/fit-count-zinb.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-count-zinb.R

zinb_full_fit.rds: R/fit-count-zinb-nuts.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-count-zinb-nuts.R

ba_gamma_fit.rds: R/fit-burn-area-gamma.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-burn-area-gamma.R

ba_lognormal_fit.rds: R/fit-burn-area-lognormal.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-burn-area-lognormal.R

ba_pareto_fit.rds: R/fit-burn-area-pareto.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-burn-area-pareto.R

ba_tpareto_fit.rds: R/fit-burn-area-tpareto.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-burn-area-tpareto.R

ba_weibull_fit.rds: R/fit-burn-area-weibull.R data/processed/stan_d.rds
		Rscript --vanilla R/fit-burn-area-weibull.R

fig/ppc-density-funs.png data/processed/burn-area-loglik.csv data/processed/area_coverage.csv: R/burn-area-model-comps.R \
	ba_gamma_fit.rds ba_lognormal_fit.rds ba_pareto_fit.rds \
	ba_tpareto_fit.rds ba_weibull_fit.rds
		Rscript --vanilla R/burn-area-model-comps.R

fig/burn-area-effs.png data/processed/burn-area-beta.csv: ba_lognormal_fit.rds R/burn-area-plots.R
		Rscript --vanilla R/burn-area-plots.R

data/processed/count-loglik.csv fig/ppc-counts.png: R/count-model-comps.R \
	pois_fit.rds nb_fit.rds zip_fit.rds zinb_fit.rds
		Rscript --vanilla R/count-model-comps.R

fig/all-coefs.png fig/count-partial-effs.png data/processed/rho_beta.csv: R/count-effplots.R \
	zinb_full_fit.rds
		Rscript --vanilla R/count-effplots.R

fig/attribution-plot.png: R/interaction-plots.R zinb_full_fit.rds test_preds.rds
		Rscript --vanilla R/interaction-plots.R

count-preds.rds fig/count-preds.png data/processed/area_df.csv data/processed/count_test_intervals.csv: R/plot-predicted-counts.R \
	zinb_full_fit.rds
		Rscript --vanilla R/plot-predicted-counts.R

fig/max-preds-l3.png test_preds.rds data/processed/predicted_totals.csv data/processed/mev_intervals.csv: count-preds.rds \
	ba_lognormal_fit.rds R/mev-plots.R data/processed/area_df.csv
		Rscript --vanilla R/mev-plots.R
		
fig/maps.png: data/processed/mtbs.rds data/processed/ecoregions.rds R/plot-study-region.R
		Rscript --vanilla R/plot-study-region.R
		
fig/spline-concept.png: R/spline-concept.R
		Rscript --vanilla R/spline-concept.R
		
push_fits: 
	aws s3 cp . s3://earthlab-mjoseph/ --recursive --exclude "*" --include "*_fit.rds"

pull_fits:
	aws s3 cp s3://earthlab-mjoseph/ . --recursive --exclude "*" --include "*_fit.rds"


clean: 
		find . -type f -name '*.png' -delete
		find . -type f -name '*.pdf' -delete
		find . -type f -name '*.rds' -delete
		find . -type f -name '*.zip' -delete
		find . -type f -name '*.tif' -delete
		rm -r data/raw/cb_2016_us_nation_20m
		rm -r data/raw/mtbs_fod_pts_data
		rm -r data/raw/us_eco_l3
		rm -r us_pbg00_2007.gdb
		rm -r data/processed/climate-data
	