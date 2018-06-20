#!/bin/bash

#This file collects all BF data from each run in model_output/
#  and then collates it in model_output/I211_J211.dat, etc
#  where each run is a column in the file

#!/bin/bash

EXEC_FILE='collect.x'
PARAM_FILE='script_parameters.dat'
NAMES_FILE='script_names.dat'

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

while read NAME; do
   for((j=0;j<$RUNS;j++))
   do
       if [ ! -d $NAME ]
       then
	   mkdir -v $NAME
       fi
       RUN=` printf "${FOLDER}run%04d/${NAME}" $j `
       NEWFILE=` printf "${NAME}/run%04d" $j `   
       egrep -v "^(#|$)" $RUN > $NEWFILE
   done
   DEST="$FOLDER$NAME"
   ./collect.x $RUNS $NAME $DEST
   rm -r $NAME
done < <(cat $NAMES_FILE)