#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include "armadillo"

int main(int argc, char* argv[]){
  std::string filename;
  int lines,lhp_samples;
  if(argc!=4){
    std::cout << "Usage: Enter also 'fn lines lhp_s' the same line.\n";
    exit(1);
  }
  else{
    filename=argv[1];
    lines=atoi(argv[2]);
    lhp_samples=atoi(argv[3]);
  }

  /*Load parameter_priors information from file
    We store twice in Names[i] to rid outselves of the first field
    in the file, "UNIFORM"
   */
  std::ifstream ifile;
  arma::mat File = arma::zeros<arma::mat>(lines,2);
  std::vector<std::string> Names(lines);
  ifile.open(filename);
  for(int i=0;i<lines;i++){
    ifile >> Names[i]; ifile >> Names[i];
    ifile >> File(i,0);
    ifile >> File(i,1);
  }
  ifile.close();
  /*Construct the generic hypercube sampling in terms of an int list, 0 to lhp_samples
    We end with a lines by lhp_s matrix
   */
  arma::mat hypercube = arma::zeros<arma::mat>(lhp_samples,lines);
  arma::vec hyperlist = arma::linspace<arma::vec>(0,lhp_samples-1,lhp_samples);
  for(int i=0;i<lines;i++){
    hyperlist=shuffle(hyperlist);
    hypercube.col(i) = hyperlist;
  }

  /*Now, construct full numerical lhp sampling using the input from file
    and the previous generic sampling matrix
  */
  float init,final,dx;
  for(int i=0;i<lines;i++){
    init = File(i,0);
    final = File(i,1);
    dx = (final-init)/(lhp_samples-1);
    hyperlist = hypercube.col(i);
    hyperlist = init + dx*hyperlist;
    hypercube.col(i) = hyperlist;
  }

  /*Construct the parameter files for each sampling 
   */
  std::ofstream ofile;
  for(int i=0;i<lhp_samples;i++){
    char fn[50];
    sprintf(fn,"../model_output/run%04d/parameters.dat",i);
    ofile.open(fn);
    for(int j=0;j<lines;j++){
      ofile << Names[j] << " " << hypercube(i,j) << '\n';
    }
    ofile.close();
  }
  
  return 0;
}
