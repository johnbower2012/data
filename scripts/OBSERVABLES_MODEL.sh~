#!/bin/bash

EXEC_FILE='observables.x'
PARAM_FILE='script_model_parameters.dat'
NAMES_FILE='script_model_names.dat'

SOURCE_FOLDER=''
DEST_FOLDER=''
LINES=0
RUNS=0
FILES=0

while read NAME VALUE; do
   if [ ${NAME} == 'SOURCE_FOLDER' ]
   then
	SOURCE_FOLDER=${VALUE}
   elif [ ${NAME} == 'DEST_FOLDER' ]
   then
	DEST_FOLDER=${VALUE}
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
   dest="$SOURCE_FOLDER$name"
   array[$i]="$dest"
done < <(egrep -v '^(#|$)' $NAMES_FILE)

echo "OBSERVABLES.sh Running:"
echo $EXEC_FILE $DEST_FOLDER $LINES $RUNS $FILES ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]}

./$EXEC_FILE $DEST_FOLDER $LINES $RUNS $FILES ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]}
