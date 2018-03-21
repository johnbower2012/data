#!/bin/bash

outfilename="outfile.dat"
path='../model_output/'
runs=`ls -d ${path}run* | wc -l`
ls ../model_output/run0000 | sort | tail -n +3 | head -6 > names.dat
names=`wc -l < names.dat`
declare -a array
for((i=1;i<$names+1;i++))
do
    name=`cut -d $'\n' -f $i names.dat`
    dest="$path$name"
    lines=`wc -l < $dest`
    array[$i]="$dest"
done

./main.x 5 5 ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} $outfilename

rm names.dat
