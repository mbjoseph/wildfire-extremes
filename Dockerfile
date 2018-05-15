FROM rocker/tidyverse:latest

MAINTAINER Max Joseph maxwell.b.joseph@colorado.edu


# rstan installation taken from jrnold/rstan
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
  	apt-utils \
  	ed \
  	libnlopt-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/

# Global site-wide config -- neeeded for building packages
RUN mkdir -p $HOME/.R/ \
  && echo "CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -flto -ffat-lto-objects  -Wno-unused-local-typedefs \n" >> $HOME/.R/Makevars

# Config for rstudio user
RUN mkdir -p $HOME/.R/ \
  && echo "CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -flto -ffat-lto-objects  -Wno-unused-local-typedefs -Wno-ignored-attributes -Wno-deprecated-declarations\n" >> $HOME/.R/Makevars \
  && echo "rstan::rstan_options(auto_write = TRUE)\n" >> /home/rstudio/.Rprofile \
  && echo "options(mc.cores = parallel::detectCores())\n" >> /home/rstudio/.Rprofile

# Install rstan
RUN install2.r --error --deps TRUE rstan




# other deps once rstan installed
RUN apt-get update && apt-get install -y --no-install-recommends \
  gdal-bin \ 
  libgdal-dev \
  libhdf4-alt-dev \
  libhdf5-dev \
  liblwgeom-dev \
  libnetcdf-dev \
  libsqlite3-dev \
  libssh2-1-dev \
  libssl-dev \
  libudunits2-dev \
  python-pip \
  texlive-full

RUN install2.r --error \
  assertthat \
  bookdown \
  cowplot \
  foreign \
  ggrepel \
  ggridges \
  ggthemes \
  hrbrthemes \
  kableExtra \
  lubridate \
  maptools \
  ncdf4 \
  pbapply \
  raster \
  RCurl \
  rgdal \
  sf \
  spdep \
  snowfall \
  viridis \
  zoo \
  ## from bioconductor
  && R -e "source('https://bioconductor.org/biocLite.R'); biocLite('BiocInstaller')" \
  && R -e "BiocInstaller::biocLite('rhdf5')" \
  ## from GitHub
  && R -e "devtools::install_github('thomasp85/patchwork')"

RUN pip install wheel \ 
  && pip install awscli
  
WORKDIR /home/rstudio/wildfire-extremes

COPY . .

RUN chown rstudio . -R

