#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include "armadillo"

int main(int argc, char* argv[]){
  std::string sourcefilename, destinationfilename;
  char tempfilename[50];
  int runs,lines;
  if(argc != 4){
    std::cout << "Usage: Enter also 'runs source_fn destination_fn'.\n";
    exit(1);
  }
  else{
    runs = atoi(argv[1]);
    sourcefilename = argv[2];
    destinationfilename = argv[3];
  }
  std::ofstream ofile;
  std::ifstream ifile;
  lines=0;

  //Determine how many lines exist
  sprintf(tempfilename, "./%s/run%04d", sourcefilename.c_str(), 0);
  ifile.open(tempfilename);
  std::string linefromfile;
  while(!ifile.eof()){
    std::getline(ifile,linefromfile);
    lines++;
  }
  ifile.close();

  //Assemble all Balance Function data
  arma::mat BF = arma::zeros<arma::mat>(lines,runs);
  arma::vec delY = arma::zeros<arma::vec>(lines);
  double number;
  for(int i=0;i<runs;i++){
    sprintf(tempfilename, "%s/run%04d", sourcefilename.c_str(), i);
    ifile.open(tempfilename);
    if(ifile){
      for(int j=0;j<lines;j++){
	ifile >> delY(j);
	ifile >> BF(j,i);
	ifile >> number;
      }
    }
    else{
      printf("%d%s",i," cannot open.");
    }
    ifile.close();
  }

  //print out the Balance Function Data
  ofile.open(destinationfilename);
  for(int i=0;i<lines;i++){
    ofile << std::setw(15) << delY(i);
    for(int j=0;j<runs;j++){
      ofile << std::setw(15) << BF(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
  
  return 0;
}
