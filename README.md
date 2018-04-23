# Spatiotemporal modeling of wildfire extremes

[![Docker Automated build](https://img.shields.io/docker/automated/mbjoseph/wildfire-extremes.svg)](https://hub.docker.com/r/mbjoseph/wildfire-extremes/)
[![Docker Build Status](https://img.shields.io/docker/build/mbjoseph/wildfire-extremes.svg)](https://hub.docker.com/r/mbjoseph/wildfire-extremes/)


This repository contains code to build spatiotemporal models of wildfire extremes in the contiguous United States. 

## Hardware requirements

We recommend at least 4 physical CPUs and 30 GB of RAM. 

## Reproducing the analysis

### Spinning up the computational environment

We have provided a Docker container that bundles up the software dependencies for this project, and provides an RStudio server instance that can be used in a web browser. 
To launch the container, run the following:

```bash
docker run -d -p 8787:8787 mbjoseph/wildfire-extremes
```

Then, navigate to port 8787 on a web browser (e.g., localhost:8787) and log in with username `rstudio`, password `rstudio`. 

### Running the analysis

To run everything, you can type the following command from a terminal: 

```
make
```

## Overview of workflow

We split the workflow into three steps: 

### 1. Data processing 

Climate, housing, and wildfire data are acquired and cleaned via:

```bash
make data
```

### 2. Model training and evaluation

TODO: describe

### 3. Manuscript generation

TODO: describe
