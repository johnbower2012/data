#include "observables.hpp"

double median_Width(arma::mat function, double factor){
  int rows = function.n_rows;
  int cols = function.n_cols;
  double median_width, sum=0.0;
  arma::vec tally = arma::zeros<arma::vec>(rows-1);
  for(int i=0;i<rows-1;i++){
    tally(i) = (function(i+1,1) + function(i,1))*(function(i+1,0) - function(i,0))/2.0;
    sum += tally(i);
    tally(i) = sum;
  }
  sum*=factor;
  for(int i=0;i<rows-1;i++){
    if(tally(i)>sum){
      median_width = (function(i+1,0) + function(i,0))/2.0;
      break;
    }
  }
  return median_width;
}

int main(int argc, char* argv[]){

  std::string *infilename,outfilename;
  int observables=7,lines,runs;
  infilename = new std::string[observables];

  if(argc<5+observables){
    std::cout << "Improper input. Enter also 'lines runs observables ifn*[7] ofn' on same line." << std::endl;
    exit(1);
  }
  else{
    infilename = new std::string[observables];
    lines=atoi(argv[1]);
    runs=atoi(argv[2]);
    observables=atoi(argv[3]);
    for(i=0;i<observables;i++){
      infilename[i]=argv[4+i];
    }
    outfilename = argv[4+observables];
  }

  
  
  return 0;
}
