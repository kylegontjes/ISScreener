#!/bin/sh

#1 = Path to results folder (can be absolute, relate, or just current folder)
#2 = Name of the final report (e.g. "2024_06_03_all_runs")

# Curate ISfinder results 
cd $1

# First report
runs=$(ls -d */ | sort | sed 's/\///' )
first_run=$(echo "$runs" | head -n1)

# 0. Create master file
master_file=$2\_ISFinder_master_report.txt
touch $master_file
# 1. Add first line of first isolate
cat $first_run/$first_run\_master.txt | head -n1 >> $master_file

#2. Add all lines from each directory's repot
for run in `echo "$runs"`
do 
echo $run
tail -n +2 $run/$run\_master.txt >> $master_file
done 