#!/bin/bash

#This file collects all BF data from each run in model_output/
#  and then collates it in model_output/I211_J211.dat, etc
#  where each run is a column in the file

path='../model_output/'
runs=`ls -d ${path}run* | wc -l`
ls ../model_output/run0000 | sort | tail -n +3 | head -6 > names.dat
lines=`wc -l < names.dat`
for((i=1;i<$lines+1;i++))
do
    names=`cut -d $'\n' -f $i names.dat`
    for((j=0;j<$runs;j++))
    do
	run=` printf "${path}run%04d/$names" $j`
	if [ ! -d $names ]
	then
	    mkdir $names
	fi
	newfile=` printf "${names}/run%04d" $j`
	cat $run | tail -n +7 > $newfile
    done
    dest="$path$names"
    ./collect.x $runs $names $dest
    rm -r $names
done

rm names.dat
