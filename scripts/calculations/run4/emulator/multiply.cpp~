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
  std::string filename1, filename2, outfile;
  if(argc<4){
    printf("Improper usage. Please enter also 'filename1 filename2 outfilename' on same line.\n");
    exit(1);
  } else{
    filename1 = argv[1];
    filename2 = argv[2];
    outfile = argv[3];
  }

  int exp_lines=1;
  int observables=12;
  int runs=1000;
  
  arma::mat exp = arma::zeros<arma::mat>(exp_lines,observables);
  arma::mat eigvecs = arma::zeros<arma::mat>(observables,observables);
  load_file(filename1, exp);
  load_file(filename2, eigvecs);
  arma::mat result = exp*eigvecs;
  write_file(outfile, result);

  return 0;
}
