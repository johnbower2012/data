/* COMMENTS SECTION:


	***********************************************************/

#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include<fstream>
#include<string>
#include "time.h"
#include "armadillo"
#include "gaussianprocess.cpp"

std::ofstream ofile;
std::ifstream ifile;

arma::mat construct_latinhypercube_sampling(int samples, arma::mat range);
void load_data_file(std::string fileName, arma::mat &X, arma::mat &Y){
  int param = X.n_cols;
  int obs = Y.n_cols;
  int train = X.n_rows;

  std::ifstream ifile;
  ifile.open(fileName);
  for(int i=0;i<train;i++){
    for(int j=0;j<param;j++){
      ifile >> X(i,j);
    }
    for(int j=0;j<obs;j++){
      ifile >> Y(i,j);
    }
  }
  ifile.close();
}
void load_beta_file(std::string betaName, arma::mat &beta){
  int load = beta.n_rows;
  int param = beta.n_cols;
  std::string temp;
  std::ifstream ifile;
  ifile.open(betaName);
  ifile >> temp;
  for(int i=0;i<load;i++){
    for(int j=0;j<param;j++){
      ifile >> beta(i,j);
    }
  }
  ifile.close();
}


int main(int argc, char* argv[]){
	//Variables for use throughout the program
  int i,j,
    train, test, 
    param, observables, hyperp,
    samples;
  unsigned seed;
  double 	x, epsilon;
  double dx, x_init, x_final;
  double sigma_n, sigma_f, l;

  //For writing to file
  std::string param_filename, obs_filename, outfilename, trainfilename;

  //Test to ensure the proper input. Terminate program if command line inputs are not given
  if(argc!=5){
    std::cout << "Improper usage. Please also enter 'filename sigma_f l sigma_n' on same line." << std::endl;
    exit(1);
  }
  else{
    param_filename = argv[1];
    sigma_f = atof(argv[2]);
    l = atof(argv[3]);
    sigma_n = atof(argv[4]);
  }
  
  /*SYSTEM and OUTPUT--
***********************************************************/
  
  train = 1000;
  test = 719;
  
  samples = 1;
  param = 4;
  observables = 12;

  epsilon = 1e-8;

  seed = std::chrono::system_clock::now().time_since_epoch().count();

  outfilename = "testset.dat";
  trainfilename = "trainset.dat";
  
  /*HYPERPARAMETERS--
***********************************************************/

  hyperp = 3;

  //hyperp_vec(0) = sigma_f;
  //hyperp_vec(1) = l;
  //hyperp_vec(2) = sigma_n;

  /*SAMPLING RANGE--
***********************************************************/
  //random number generator
  std::default_random_engine generator(seed);
  std::normal_distribution<double> dist(0,1.0);

  //armadillo matrices and vectors for computational use
  arma::mat xvec_train_mat = arma::zeros<arma::mat>(train,param),
    yvec_train_mat = arma::zeros<arma::mat>(train,observables),
    xvec_test_mat = arma::zeros<arma::mat>(test,param),
    
    output_mat = arma::zeros<arma::mat>(test,param + 3 + samples);
  
  arma::vec hyperp_vec = arma::zeros<arma::vec> (hyperp);
  
  //hyperparamter vector for use in 'function'
  hyperp_vec(0) = sigma_f;
  hyperp_vec(1) = l;
  hyperp_vec(2) = sigma_n;

  load_data_file(param_filename,xvec_train_mat,yvec_train_mat);
  std::string betaName = "beta.dat";
  arma::mat beta = arma::zeros<arma::mat>(param*observables,2);
  load_beta_file(betaName,beta);

  observables=1;
  yvec_train_mat.resize(train,1);

  float init,final;
  arma::mat range = arma::zeros<arma::mat>(param,2);  
  for(int i=0;i<param;i++){
    init = range(i,0) = 0.01;
    final = range(i,1) = 2.0;
  }
  xvec_test_mat = construct_latinhypercube_sampling(test, range);


/***********************************************************/
  /*END SETUP--*/
    /*BEGIN COMPUTATION
***********************************************************/
  
  /////////////////////////
  /*

  train = 17;
  test = 23;
  param = 1;
  observables=1;
  samples=1;
  xvec_train_mat = arma::zeros<arma::mat>(train,1);
  yvec_train_mat = arma::zeros<arma::mat>(train,1);
  xvec_test_mat = arma::zeros<arma::mat>(test,1);
  for(int i=0;i<train;i++){
    xvec_train_mat(i,0) = ((double) (i+1))/train*(10.0);;
    //yvec_train_mat(i,0) = xvec_train_mat(i,0);
    yvec_train_mat(i,0) = sin(xvec_train_mat(i,0));
  }
  for(int i=0;i<test;i++){
    xvec_test_mat(i,0) = ((double) (i+1))/test*(14.0) - 2.0;
  }
  */
  ///////////////////////

  //Execute posterior function generation
  //  output_mat = gaussian_process_solver_basic(kernel_square_exp_function, xvec_train_mat, yvec_train_mat, xvec_test_mat, samples, hyperp_vec, epsilon);
  output_mat = gaussian_process_solver_regression(kernel_square_exp_function, xvec_train_mat, yvec_train_mat, xvec_test_mat, samples, hyperp_vec, epsilon, regression_linear_function,beta);
  //Write the output file
  write_output(output_mat, test, param, observables, samples, outfilename);
  //Write the trainingset to file
  write_trainset(xvec_train_mat, yvec_train_mat, trainfilename);

  return 0;

}

arma::mat construct_latinhypercube_sampling(int samples, arma::mat range){
  int parameters = range.n_rows;
  arma::mat hypercube = arma::zeros<arma::mat>(samples,parameters);
  arma::vec hyperlist = arma::linspace<arma::vec>(0,samples-1,samples);
  for(int i=0;i<parameters;i++){
    hyperlist=shuffle(hyperlist);
    hypercube.col(i) = hyperlist;
  }

  float init,final,dx;
  for(int i=0;i<parameters;i++){
    init = range(i,0);
    final = range(i,1);
    dx = (final-init)/(samples-1);
    hyperlist = hypercube.col(i);
    hyperlist = init + dx*hyperlist;
    hypercube.col(i) = hyperlist;
  }
  return hypercube;
}
