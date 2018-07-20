#include<iostream>
#include<iomanip>
#include<cmath>
#include "emulator_class.cpp"
#include<armadillo>

int main(int argc, char* argv[]){
  int train = 1000;
  int test = 719;
  int param = 4;
  int observables=12;
  int hp = 3;
  double epsilon=1e-8;

  std::string fileName, betaName,hypName;

  arma::mat X = arma::zeros<arma::mat>(train,param);
  arma::mat X_s = arma::zeros<arma::mat>(test,param);
  arma::mat range = arma::zeros<arma::mat>(param,2);
  arma::mat H = arma::zeros<arma::mat>(observables,hp);
  arma::mat beta = arma::zeros<arma::mat>(observables,param+1);
  arma::mat Y = arma::zeros<arma::mat>(train,observables);
  arma::mat y = arma::zeros<arma::mat>(train,1);

  if(argc<4){
    printf("Improper usage. Please enter 'fileName betaName hyperparametername' on same line.\n");
    exit(1);
  } else{
    fileName = argv[1];
    betaName = argv[2];
    hypName = argv[3];
  }
  
  load_data_file(fileName,X,Y);
  load_beta_file(betaName,beta);
  load_file(hypName,H);

  for(int i=0;i<param;i++){
    range(i,0) = 0.01;
    range(i,1) = 2.0;
  }
  X_s = construct_latinhypercube_sampling(test,range);

  emulator gauss(X,H,beta,epsilon);
  arma::mat output = gauss.emulate(X_s,Y);

  write_output(output,param,observables,"test.dat");
  write_trainset(X, Y, "train.dat");
  return 0;
}
