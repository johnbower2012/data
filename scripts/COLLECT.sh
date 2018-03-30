#!/bin/bash

#This file collects all BF data from each run in model_output/
#  and then collates it in model_output/I211_J211.dat, etc
#  where each run is a column in the file

path='../model_output/'
runs=`ls -d ${path}run* | wc -l`
find ../model_output/run0000/I* | cut -d '/' -f 4 > names.dat
printf "allcharges.dat\n" >> names.dat
names=`wc -l < names.dat`
for((i=1;i<$names+1;i++))
do
    name=`cut -d $'\n' -f $i names.dat`
    for((j=0;j<$runs;j++))
    do
	run=` printf "${path}run%04d/$name" $j`
	if [ ! -d $name ]
	then
	    mkdir -v $name
	fi
	newfile=` printf "${name}/run%04d" $j `
	grep -v "^#" $run > $newfile
    done
    dest="$path$name"
    echo "$dest finished"
    ./collect.x $runs $name $dest
    rm -r $name
done

rm names.dat
