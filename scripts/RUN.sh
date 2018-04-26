#!/bin/bash

for((i=$1;i<$2;i++))
do
    fn=` printf "run%04d" $i`
    echo -e "\n+++++++++ Running Balance on $fn +++++++++\n"
    ./balance $fn
done

echo Finished...\n
