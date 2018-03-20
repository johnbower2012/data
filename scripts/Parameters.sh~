#!/bin/bash

#This script will create all necessary run folders and
#  parameter files required to run 'balance'

#File location
file_priors='../model_output/parameter_priors.dat'
file_fixed='../model_output/fixed_parameters.dat'

#How many parameters will we conduct the sampling?
lines=`cat $file_priors | cut -d ' ' -f 3,4 | wc -l`

#How many samples? First input is samples. Default 1000.
if [ $# -lt 1 ]
then
    latinHyperCube_samples=1000
else
    latinHyperCube_samples=$1
fi

for((i=0;i<latinHyperCube_samples;i++))
    do
	fn=`printf "../model_output/run%04d" $i`
	if [ ! -d $fn ]
	then
	    mkdir -pv $fn
	fi
	fn=`printf "../model_output/run%04d/fixed_parameters.dat" $i`
	if [ ! -d $fn ]
	then
	    cp $file_fixed $fn
	fi
    done

#Execute parameters.x
./parameters.x $file_priors $lines $latinHyperCube_samples
