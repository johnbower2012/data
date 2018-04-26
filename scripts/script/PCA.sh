#!/bin/bash

EXEC_FILE='pca.x'
OUTFILENAME='pca.dat'

PARAM_FILE='script_parameters.dat'
NAMES_FILE='script_names.dat'

FOLDER=''
LINES=0
RUNS=0
FILES=0
while read NAME VALUE; do
   if [ ${NAME} == 'FOLDER' ]
   then
	FOLDER=${VALUE}
   elif [ ${NAME} == 'LINES' ]
   then
       LINES=${VALUE}
   elif [ ${NAME} == 'RUNS' ]
   then 
       RUNS=${VALUE}
   elif [ ${NAME} == 'FILES' ]
   then
       FILES=${VALUE}
   else
       echo "Variable unknown: $NAME $VALUE"
   fi
done < <(egrep -v '^(#|$)' $PARAM_FILE)

declare -a array
for((i=1;i<$FILES+1;i++))
do
    NAME=` cut -d $'\n' -f $i $NAMES_FILE `
    DEST=$FOLDER$NAME
    array[$i]=$DEST
done

./pca.x $LINES $RUNS $FILES ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]} $OUTFILENAME
