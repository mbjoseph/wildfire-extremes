#!bin/bash
# this script splits the big precip data into small year-specific csv files

for i in `seq 1991 2013`;
do
	grep $i data/raw/us-precip.csv > data/processed/precip-$i.csv
done

