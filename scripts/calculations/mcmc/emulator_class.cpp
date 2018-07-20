#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include "time.h"
#include "armadillo"

class emulator{
public:
  int train_points;
  int parameter_count;
  int observable_count;
  int hyperparameter_count;

  double epsilon;
  arma::vec noise;

  arma::mat X;
  arma::mat *kernel;
  arma::mat *kernel_inverse;
  arma::mat hyperparameters;
  arma::mat *H;
  arma::mat beta;

  emulator(arma::mat X, arma::mat hyperparameters, arma::mat beta, double epsilon);

  arma::mat kernel_function(arma::mat A, arma::mat B, int obs_index);
  arma::mat regression_linear_function(arma::mat X_mat, int obs_index);
  arma::mat emulate(arma::mat X_s, arma::mat Y_mat);
  arma::mat emulate_nr(arma::mat X_s, arma::vec y);
};



emulator::emulator(arma::mat new_X, arma::mat new_hyperparameters, arma::mat new_beta, double new_epsilon){
  this->X = new_X;
  this->hyperparameters = new_hyperparameters;
  this->beta = new_beta;
  this->epsilon = new_epsilon;

  this->train_points = X.n_rows;
  this->parameter_count = X.n_cols;
  this->observable_count = new_beta.n_rows;
  this->hyperparameter_count = hyperparameters.n_cols;

  this->noise = hyperparameters.col(hyperparameter_count - 1);

  arma::mat I_train = arma::eye<arma::mat>(train_points,train_points);
  kernel = new arma::mat[observable_count];
  kernel_inverse = new arma::mat[observable_count];
  H = new arma::mat[observable_count];
  for(int i=0;i<this->observable_count;i++){
    this->kernel[i] = kernel_function(X,X,i) + I_train*(noise(i)*noise(i)+epsilon);
    this->kernel_inverse[i] = kernel[i].i();
    this->H[i] = regression_linear_function(X,i);
  }
}
arma::mat emulator::kernel_function(arma::mat A, arma::mat B, int obs_index){
  double xi, xj, sum, diff,
    theta1, theta2;
  int i,j,k,
    ai = A.n_rows, bi = B.n_rows;
  arma::mat kernel_matrix = arma::zeros<arma::mat>(ai,bi);

  theta1 = hyperparameters(obs_index,0)*hyperparameters(obs_index,0);
  theta2 = 1.0/hyperparameters(obs_index,1)/hyperparameters(obs_index,1)/2.0;

  for(i=0;i<ai;i++){
    for(j=0;j<bi;j++){
      sum = 0.0;
      for(k=0;k<parameter_count;k++){
	xi = A(i,k);
	xj = B(j,k);
	diff = xi - xj;
	sum += diff*diff;
      }
      kernel_matrix(i,j) = theta1*exp(-sum*theta2);
    }
  }

  return kernel_matrix;
}
arma::mat emulator::regression_linear_function(arma::mat X_mat, int obs_index){
  int points = X_mat.n_rows;
  arma::mat H_mat = arma::zeros<arma::mat>(1,points);

  for(int i=0;i<points;i++){
    H_mat(0,i) = beta(obs_index,0);
    for(int j=0;j<parameter_count;j++){
      H_mat(0,i) += beta(obs_index,j+1)*X_mat(i,j);
    }
  }
  return H_mat;
}
arma::mat emulator::emulate(arma::mat X_s, arma::mat Y){
  //random number generator
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> dist(0,1.0);

  int i, j,
    test_points = X_s.n_rows;

  //armadillo matrices and vectors for computational use
  arma::mat kernel_s = arma::zeros<arma::mat>(train_points,test_points),
    kernel_ss = arma::zeros<arma::mat>(test_points,test_points),
    
    mean_mat = arma::zeros<arma::mat>(test_points,1),
    var_mat = kernel_ss,
    L = var_mat,

    H_s = arma::zeros<arma::mat>(1,test_points),
    R = arma::zeros<arma::mat>(train_points,test_points),
    beta_mat = arma::zeros<arma::mat>(train_points,1),

    temp_pp = arma::zeros<arma::mat>(1,1),
    temp_pp_inverse = arma::zeros<arma::mat>(1,1),
				
    post_func = arma::zeros<arma::mat>(test_points,1),
    random_sample = post_func,

    I_test = arma::eye<arma::mat>(test_points,test_points),

    output_mat = arma::zeros<arma::mat>(test_points, parameter_count + observable_count*3);


  for(int i=0;i<observable_count;i++){
    //calculate kernel matrices
    kernel_s = kernel_function(X,X_s,i);
    kernel_ss = kernel_function(X_s,X_s,i) + I_test*epsilon;

    //numerical stability factor in solving the Cholesky Decomposition
    kernel_ss += I_test*epsilon;

    //calculate regression related matrices, R and beta
    H_s = regression_linear_function(X_s,i);

    R = H_s - H[i]*kernel_inverse[i]*kernel_s;

    temp_pp = H[i]*kernel_inverse[i]*H[i].t();
    temp_pp_inverse = temp_pp.i();
    beta_mat = temp_pp_inverse*H[i]*kernel_inverse[i]*Y.col(i);
  
    //calculate MEAN of distribution
    mean_mat = kernel_s.t()*kernel_inverse[i]*Y.col(i) + R.t()*beta_mat;
 
    //Calculate COVARIANCE of distribution
    var_mat = kernel_ss - kernel_s.t()*kernel_inverse[i]*kernel_s + R.t()*temp_pp_inverse*R;	

    //Cholesky Decomposition, 'sqrt' of Variance matrix, here the lower triangular form
    L = arma::chol(var_mat,"lower");

    //generator random distribution
    //create matrix of "Mean" values for easy generation of output function values

    for(j=0;j<test_points;j++){
      random_sample(j,0) = dist(generator);
    }

    //sample posterior functions
    post_func = mean_mat + L*random_sample;

    output_mat.col(parameter_count + i) = mean_mat.col(0);
    output_mat.col(parameter_count + observable_count + i) = var_mat.diag();
    output_mat.col(parameter_count + 2*observable_count + i) = post_func.col(0);
  }

  //generate output matrix:
  //		param1, param2, ... , mean1, var1, mean2, var2, ... , post_func1, post_func2, ...
  for(i=0;i<test_points;i++){
    for(j=0;j<parameter_count;j++){
      output_mat(i,j) = X_s(i,j);
    }
  }
  
  return output_mat;

}
/*
arma::mat emulator::emulate_nr(arma::mat X_s, arma::vec y){
  //random number generator
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> dist(0,1.0);

  int i, j,
    test_points = X_s.n_rows;
  
  double noise = hyperparameters(hyperparameter_count-1)*hyperparameters(hyperparameter_count-1);

  //armadillo matrices and vectors for computational use
  arma::mat kernel_s = arma::zeros<arma::mat>(train_points,test_points),
    kernel_ss = arma::zeros<arma::mat>(test_points,test_points),

    mean_mat = arma::zeros<arma::mat>(test_points,1),
    var_mat = kernel_ss,
    L = var_mat,
				
    post_func = arma::zeros<arma::mat>(test_points,1),
    random_sample = post_func,

    I_train = arma::eye<arma::mat>(train_points,train_points),
    I_test = arma::eye<arma::mat>(test_points,test_points),

    output_mat = arma::zeros<arma::mat>(test_points, parameter_count + 2 + 1);



  //calculate kernel matrices
  kernel_s = kernel_function(X,X_s);
  kernel_ss = kernel_function(X_s,X_s);

  //numerical stability factor in solving the Cholesky Decomposition
  kernel_ss += I_test*epsilon;

  //calculate MEAN of distribution
  mean_mat = kernel_s.t()*kernel_inverse*y;
  //Calculate COVARIANCE of distribution
  var_mat = kernel_ss - kernel_s.t()*kernel_inverse*kernel_s;
  //Cholesky Decomposition, 'sqrt' of Variance matrix, here the lower triangular form
  L = arma::chol(var_mat,"lower");

  //generator random distribution
  //create matrix of "Mean" values for easy generation of output function values
  for(j=0;j<test_points;j++){
    random_sample(j,0) = dist(generator);
  }
 
  //sample posterior functions
  post_func = mean_mat + L*random_sample;

  //generate output matrix:
  //		param1, param2, ... , mean1, var1, mean2, var2, ... , post_func1, post_func2, ...
  for(i=0;i<test_points;i++){
    for(j=0;j<parameter_count;j++){
      output_mat(i,j) = X_s(i,j);
    }
    output_mat(i,parameter_count) = mean_mat(i,0);
    output_mat(i,parameter_count+1) = var_mat(i,0);
    output_mat(i,parameter_count+2) = post_func(i,0);
  }

  return output_mat;
}
*/
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

void write_output(arma::mat output_mat, int param, int obs, std::string outfilename){
  //write to file as:
  //		param1, param2, ... , mean value, mean-2*std, mean+2*std, output1, output2,....
  int i, j, 
    test_points = output_mat.n_rows;
  std::ofstream ofile;
  
  ofile.open(outfilename);
  for(i=0;i<test_points;i++){
    for(j=0;j<param;j++){
      ofile << " " << output_mat(i,j);
    }
    for(j=0;j<obs;j++){
      ofile << " " << output_mat(i,param+j);
      ofile << " " << output_mat(i,param+j) + 2.0*sqrt(output_mat(i,param+obs+j));
      ofile << " " << output_mat(i,param+j) - 2.0*sqrt(output_mat(i,param+obs+j));
      ofile << " " << output_mat(i,param+2*obs+j);
    }
    ofile << std::endl;
  }
  ofile.close();
}
void write_trainset(arma::mat X_mat, arma::mat Y_mat, std::string outfilename){
  //Write to file as:
  //		param1, param2, ... , training value
  int i, j,
    train = X_mat.n_rows,
    observables = Y_mat.n_cols,
    param = X_mat.n_cols;
  std::ofstream ofile;

  ofile.open(outfilename);
  for(i=0;i<train;i++){
    for(j=0;j<param;j++){
      ofile << " " << X_mat(i,j);
    }
    for(j=0;j<observables;j++){
      ofile << " " << Y_mat(i,j);
    }
    ofile << std::endl;
  }
  ofile.close();
}
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
void write_csv_file(std::string fileName, std::string title, arma::mat &file){
  int write = file.n_rows;
  int param = file.n_cols;
  std::string temp;
  std::ofstream ofile;
  ofile.open(fileName);
  ofile << title << std::endl;
  for(int i=0;i<write;i++){
    for(int j=0;j<param;j++){
      ofile << " " <<  file(i,j);
    }
    ofile << std::endl;
  }
  ofile.close();
}
