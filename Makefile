data-dir = data/processed
gdb = data/raw/us_pbg00_2007.gdb
figs = fig/figure_7.pdf \
	fig/figure_8.pdf fig/figure_4.pdf \
	fig/figure_6.pdf \
	fig/figure_5.pdf fig/figure_11.pdf \
	fig/figure_9.pdf \
	fig/figure_2.pdf \
	fig/figure_3.pdf \
	fig/figure_10.pdf \
	fig/figure_12.pdf
	
tables = data/processed/burn-area-loglik.csv data/processed/count-loglik.csv \
	data/processed/rho_beta.csv data/processed/count_test_intervals.csv \
	data/processed/area_df.csv data/processed/predicted_totals.csv \
	data/processed/area_coverage.csv data/processed/mev_intervals.csv \
	data/processed/burn-area-beta.csv \
	data/processed/million-er-mon.csv

all: main.pdf appendix-s1.pdf appendix-s2.pdf tex_src.zip

main.pdf: $(figs) $(tables) main.Rmd library.bib header.sty
		Rscript -e "rmarkdown::render('main.Rmd')"
		
appendix-s1.pdf: appendix-s1.Rmd header.sty
		Rscript -e "rmarkdown::render('appendix-s1.Rmd')"
		
appendix-s2.pdf: appendix-s2.Rmd header.sty
		Rscript -e "rmarkdown::render('appendix-s2.Rmd')"

data/processed/stan_d.rds data/processed/mtbs.rds data/processed/ecoregions.rds: R/make-stan-d.R \
	$(data-dir)/ecoregion_summaries.csv \
	$(data-dir)/housing_density.csv \
	data/raw/mtbs_fod_pts_data/mtbs_fod_pts_data/mtbs_fod_pts_DD.shp
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
		wget -nc -np https://s3-us-west-2.amazonaws.com/earthlab-mjoseph/us_pbg00_2007.gdb.zip
		unzip us_pbg00_2007.gdb.zip -d data/raw/
		rm us_pbg00_2007.gdb.zip
		
data/raw/us_eco_l3/us_eco_l3.shp: 
		wget -q -nc -O data/raw/us_eco_l3.zip ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip
		unzip -o data/raw/us_eco_l3.zip -d data/raw/us_eco_l3
		rm data/raw/us_eco_l3.zip

data/raw/mtbs_fod_pts_data/mtbs_fod_pts_data/mtbs_fod_pts_DD.shp:
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

fig/figure_7.pdf data/processed/burn-area-loglik.csv data/processed/area_coverage.csv: R/burn-area-model-comps.R \
	ba_gamma_fit.rds ba_lognormal_fit.rds ba_pareto_fit.rds \
	ba_tpareto_fit.rds ba_weibull_fit.rds
		Rscript --vanilla R/burn-area-model-comps.R

fig/figure_8.pdf data/processed/burn-area-beta.csv: ba_lognormal_fit.rds R/burn-area-plots.R
		Rscript --vanilla R/burn-area-plots.R

data/processed/count-loglik.csv fig/figure_4.pdf: R/count-model-comps.R \
	pois_fit.rds nb_fit.rds zip_fit.rds zinb_fit.rds
		Rscript --vanilla R/count-model-comps.R

fig/figure_6.pdf fig/figure_5.pdf data/processed/rho_beta.csv: R/count-effplots.R \
	zinb_full_fit.rds
		Rscript --vanilla R/count-effplots.R

fig/figure_11.pdf: R/interaction-plots.R zinb_full_fit.rds
		Rscript --vanilla R/interaction-plots.R

count-preds.rds data/processed/area_df.csv data/processed/count_test_intervals.csv: R/plot-predicted-counts.R \
	zinb_full_fit.rds
		Rscript --vanilla R/plot-predicted-counts.R

fig/figure_9.pdf fig/figure_10.pdf test_preds.rds data/processed/predicted_totals.csv data/processed/million-er-mon.csv data/processed/mev_intervals.csv: count-preds.rds \
	ba_lognormal_fit.rds R/mev-plots.R
		Rscript --vanilla R/mev-plots.R
		
fig/figure_2.pdf: data/processed/mtbs.rds data/processed/ecoregions.rds R/plot-study-region.R
		Rscript --vanilla R/plot-study-region.R
		
fig/figure_3.pdf: R/spline-concept.R
		Rscript --vanilla R/spline-concept.R

fig/figure_12.pdf: R/wallow-case-study.R data/processed/ecoregions.rds data/raw/climate-data.csv
		Rscript --vanilla R/wallow-case-study.R
		
push_fits: 
	aws s3 cp . s3://earthlab-mjoseph/ --recursive --exclude "*" --include "*_fit.rds"

pull_fits:
	aws s3 cp s3://earthlab-mjoseph/ . --recursive --exclude "*" --include "*_fit.rds"

tex_src.zip: $(figs) main.pdf
	zip tex_src.zip fig/*.pdf main.tex appendix-s1.tex appendix-s2.tex


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
