#This file will run all necessary programs in a continuous fashion.

#Parameters.sh will collect the necessary parameters files and create the run directories
./PARAMETERS.sh

#Run.sh runs the program
cd ..
./RUN.sh 0 1000
cd scripts/

#Collect.sh arranges the data into a usable format
./COLLECT.sh

#PCA.sh performs the analysis
./PCA.sh
