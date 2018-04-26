#!/bin/bash

EXEC_FILE='observables.x'

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
i=0
while read name; do
   i=$(($i + 1))
   dest="$FOLDER$name"
   array[$i]="$dest"
done < <(cat $NAMES_FILE)

echo "Running:"
echo $EXEC_FILE $LINES $RUNS $FILES ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]}
echo ""

./$EXEC_FILE $LINES $RUNS $FILES ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]}