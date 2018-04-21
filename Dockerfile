FROM earthlab/mstm-aws

MAINTAINER Max Joseph maxwell.b.joseph@colorado.edu

RUN apt-get update && apt-get install -y --no-install-recommends gdal-bin

COPY . /home/wildfire-extremes/

WORKDIR /home/wildfire-extremes
