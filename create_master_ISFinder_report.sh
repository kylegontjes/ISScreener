#!/bin/sh

#1 = Path to results folder 

# Curate ISfinder results 
cd $1

## 0. Create file
filename=`basename $1`
master_file=$filename\_master.txt

touch $master_file
# 1. Add first line of first isolate
first_isolate=$(ls -d */ | sort | sed 's/\///' | head -n1)

cat $first_isolate/panISa/$first_isolate\_ISFinder.txt | head -n1 >> $master_file

# 2. Add all isolates
for isolate in `ls -d */ | sort | sed 's/\///'`
do 
echo $isolate
tail -n +2 $isolate/panISa/$isolate\_ISFinder.txt >> $master_file
done
