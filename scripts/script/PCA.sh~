#!/bin/bash

outfilename="outfile.dat"
path='../model_output_false/'
runs=`ls -d ${path}run* | wc -l`
find ${path}run0000/I* | cut -d '/' -f 4 > names.dat
printf "allcharges.dat\n" >> names.dat
names=`wc -l < names.dat`
declare -a array
for((i=1;i<$names+1;i++))
do
    name=` cut -d $'\n' -f $i names.dat `
    dest="$path$name"
    array[$i]="$dest"
done

lines=` wc -l < $dest `

./pca.x $lines $runs $names ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]} $outfilename

rm names.dat
