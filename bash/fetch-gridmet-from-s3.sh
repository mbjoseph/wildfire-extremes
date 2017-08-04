#!/bin/sh

aws s3 cp s3://earthlab-gridmet/pet ../data/processed/ --recursive
aws s3 cp s3://earthlab-gridmet/pr ../data/processed/ --recursive
aws s3 cp s3://earthlab-gridmet/vs ../data/processed/ --recursive
