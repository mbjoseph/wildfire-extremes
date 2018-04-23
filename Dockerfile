FROM earthlab/mstm-aws

MAINTAINER Max Joseph maxwell.b.joseph@colorado.edu

RUN apt-get update && apt-get install -y --no-install-recommends gdal-bin

WORKDIR /home/rstudio/wildfire-extremes

COPY . .

RUN chown rstudio . -R

