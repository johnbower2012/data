#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<armadillo>

void load_file(std::string fileName, arma::mat &file){
  int load = file.n_rows;
  int param = file.n_cols;
  std::ifstream ifile;
  ifile.open(fileName);
  for(int i=0;i<load;i++){
    for(int j=0;j<param;j++){
      ifile >> file(i,j);
    }
  }
  ifile.close();
}
void load_file(std::string fileName, arma::vec &file){
  int load = file.n_elem;
  std::ifstream ifile;
  ifile.open(fileName);
  for(int i=0;i<load;i++){
    ifile >> file(i);
  }
  ifile.close();
}
void write_file(std::string fileName, arma::mat &file){
  int write = file.n_rows;
  int param = file.n_cols;
  std::string temp;
  std::ofstream ofile;
  ofile.open(fileName);
  for(int i=0;i<write;i++){
    for(int j=0;j<param;j++){
      ofile << " " <<  file(i,j);
    }
    ofile << std::endl;
  }
  ofile.close();
}

int main(int argc, char* argv[]){
  std::string expname, errorname, meanname, eigvecname, outfilename;
  if(argc<6){
    printf("Improper usage. Please enter also 'expname errorname meanname eigvecname outfilename' on same line.\n");
    exit(1);
  } else{
    expname = argv[1];
    errorname = argv[2];
    meanname = argv[3];
    eigvecname = argv[4];
    outfilename = argv[5];
  }

  int exp_lines=1;
  int observables=12;
  int runs=1000;
  
  arma::mat exp = arma::zeros<arma::mat>(exp_lines,observables);
  arma::mat zmat = arma::zeros<arma::mat>(exp_lines,observables);
  arma::vec error = arma::zeros<arma::vec>(observables);
  arma::vec mean = arma::zeros<arma::vec>(observables);
  arma::mat eigvecs = arma::zeros<arma::mat>(observables,observables);
  load_file(expname, exp);
  load_file(errorname, error);
  load_file(meanname, mean);
  load_file(eigvecname, eigvecs);
  for(int i=0;i<exp_lines;i++){
    for(int j=0;j<observables;j++){
      zmat(i,j) = (exp(i,j) - mean(j))/error(j);
    }
  }
  arma::mat result = zmat*eigvecs;
  write_file(outfilename, result);

  return 0;
}
