#include<iostream>
#include<fstream>
#include<cmath>
#include "armadillo"

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

  std::string *infilename,*outfilename;
  int observables=7,lines,runs;
  int i,j,k;
  infilename = new std::string[observables];
  outfilename = new std::string[observables];

  if(argc<4+observables){
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
  }

  std::ifstream ifile;
  std::ofstream ofile;
  arma::mat *val_matrix;
  arma::vec delY, *obs_vec;
  val_matrix = new arma::mat[observables];
  obs_vec = new arma::vec[observables];
  delY = arma::zeros<arma::vec>(lines);
  for(i=0;i<observables;i++){
    val_matrix[i] = arma::zeros<arma::mat>(runs,lines);
    obs_vec[i] = arma::zeros<arma::vec>(2);
    ifile.open(infilename[i]);
    printf("+++++ %s +++++\n", infilename[i].c_str());
    for(j=0;j<lines;j++){
      for(k=0;k<runs;k++){
	if(k==0){
	  ifile >> delY(j);
	} else{
	  ifile >> val_matrix[i](k-1,j);
	}
      }
    }
    ifile.close();
  }

  
  /*
    for(i=0;i<observables;i++){
    delimiter="/";
    position=infilename[i].rfind(delimiter);
    if(position!=std::string::npos){
      infilename[i].erase(infilename[i].begin(),infilename[i].begin()+position+delimiter.length());
    }
    for(j=0;j<lines;j++){
      for(k=0;k<lines;k++){
	print_matrix(j,k+1) = eigvec_matrix[i](j,k);
      }
      print_matrix(j,0) = delY(j);
    }					    
    writeFile(infilename[i].c_str(),print_matrix);

    delimiter = ".";
    position=infilename[i].rfind(delimiter);
    if(position!=std::string::npos){
      infilename[i].erase(infilename[i].begin()+position,infilename[i].end());
    }
    infilename[i] += "_lambda.dat";
    writeFile(infilename[i].c_str(),eigval_vec[i]);
  }
    */

  delete[] val_matrix;
  delete[] obs_vec;
  delete[] infilename;
  delete[] outfilename;
  
  return 0;
}
