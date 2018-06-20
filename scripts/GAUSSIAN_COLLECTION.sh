#!/bin/bash

#This file collects all BF data from each run in model_output/
#  and then collates it in model_output/I211_J211.dat, etc
#  where each run is a column in the file

#!/bin/bash

COLLECT_FILE='parameters.dat'
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

NEWFILE=` printf "${DEST_FOLDER}${COLLECT_FILE}" `
if [ ! -d $DEST_FOLDER ]
then
    mkdir -v $DEST_FOLDER
fi
if [ -f $NEWFILE ]
then
    rm $NEWFILE
fi
touch $NEWFILE
echo $NEWFILE
for((j=0;j<$RUNS;j++))
do
    RUN=` printf "${SOURCE_FOLDER}run%04d/${COLLECT_FILE}" $j `
    egrep -v "^(#|$)" $RUN | tr -s ' ' | cut -d ' ' -f 2 | tr '\n' ' ' >> $NEWFILE
    printf "\n" >> $NEWFILE
done