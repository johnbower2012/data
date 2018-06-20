#include<iostream>
#include<cmath>
#include<iomanip>
#include "armadillo"

/******
       Functions to create the y~ matrix and conduct PCA
 ******/

//calculate y_tilde, with zero mean and std divided by error, i.e. [z = (y-ybar)/error]
void tilde_function(const arma::mat &matrix, const arma::vec &error, arma::vec &mean, arma::mat &tilde);
//calculate y_tilde, with mean zero and unit std, i.e. [z = (y-ybar)/ystd]
void tilde_function(const arma::mat &matrix, arma::mat &covariance, arma::vec &mean, arma::mat &tilde);
//calculate covariance of y~ matrix
void covariance_function(const arma::mat &matrix, arma::mat &covariance);
//sort eigenvectors by decreasing order of eigenvalues
void sort_eigen_function(arma::vec &eigval, arma::mat &eigvec);
// END: Functions to create the y~ matrix and conduct PCA

/******
       Functions for statistics
 ******/

double median_width_function(const arma::mat &function, double factor);
double median_width_fabs_function(const arma::mat &function, double factor);
double zeroth_moment(const arma::mat &function);
double first_moment(const arma::mat &function);
double second_moment(const arma::mat &function);

// END: Functions for statistics


int main(int argc, char* argv[]){
  arma::mat matrix = arma::zeros<arma::mat>(10,5);
  arma::mat tilde = arma::zeros<arma::mat>(10,5);
  arma::mat covariance = arma::zeros<arma::mat>(5,5);
  arma::mat eigvec_matrix = arma::zeros<arma::mat>(5,5);
  arma::vec eigval_vector = arma::zeros<arma::vec>(5);
  arma::vec mean = arma::zeros<arma::vec>(5);
  arma::vec error = arma::zeros<arma::vec>(5);

  for(int i=0;i<5;i++){
    for(int j=0;j<10;j++){
      matrix(j,i) = (i+1)*(j+1)*(j+1);
    }
    error(i) = 0.10;
  }
  arma::mat function(1000,2);
  for(int i=0;i<1000;i++){
    function(i,0) = (i+1)/100.0;
    function(i,1) = (i+1)/100.0;
  }
  printf("%f\n",median_width_function(function,0.5));
  printf("%f\n",median_width_function(function,0.5));
  printf("%f\n",zeroth_moment(function));
  printf("%f\n",first_moment(function));
  printf("%f\n",second_moment(function));
  printf("%f\n",10.0/sqrt(2));
  
  tilde_function(matrix,error,mean,tilde);
  covariance_function(tilde,covariance);
  arma::eig_sym(eigval_vector, eigvec_matrix, covariance);
  sort_eigen_function(eigval_vector,eigvec_matrix);

  std::cout << "matrix:\n";
  matrix.print();
  std::cout << "eigval_vector:\n";
  eigval_vector.print();
  std::cout << "eigvec_matrix:\n";
  eigvec_matrix.print();

 
  tilde_function(matrix,covariance,mean,tilde);
  covariance_function(tilde,covariance);
  arma::eig_sym(eigval_vector, eigvec_matrix, covariance);
  sort_eigen_function(eigval_vector,eigvec_matrix);
  
  std::cout << "matrix:\n";
  matrix.print();
  std::cout << "eigval_vector:\n";
  eigval_vector.print();
  std::cout << "eigvec_matrix:\n";
  eigvec_matrix.print();
  
  return 0;
}

/******

       Functions to create the y~ matrix and conduct PCA

 ******/

//calculate y_tilde, with zero mean and std divided by error, i.e. [z = (y-ybar)/error]
void tilde_function(const arma::mat &matrix, const arma::vec &error, arma::vec &mean, arma::mat &tilde){
  int repetitions = matrix.n_rows;
  int observables = matrix.n_cols;
  mean = arma::zeros<arma::vec>(observables);
  tilde = arma::zeros<arma::mat>(repetitions,observables);

  for(int i=0;i<observables;i++){
    for(int j=0;j<repetitions;j++){
      mean(i) += matrix(j,i);
    }
    mean(i) /= (double) repetitions;
  }
  for(int i=0;i<repetitions;i++){
    for(int j=0;j<observables;j++){
      tilde(i,j) = (matrix(i,j) - mean(j))/(error(j)*matrix(i,j));
    }
  }
}
//calculate y_tilde, with mean zero and unit std, i.e. [z = (y-ybar)/ystd]
void tilde_function(const arma::mat &matrix, arma::mat &covariance, arma::vec &mean, arma::mat &tilde){
  int i, j, k,
    repetitions = matrix.n_rows,
    observables = matrix.n_cols;
  double xi_mean, xj_mean, xi_std;
  
  tilde = arma::zeros<arma::mat> (repetitions,observables);
  covariance = arma::zeros<arma::mat> (observables,observables);
  mean = arma::zeros<arma::vec>(observables);
  for(i=0;i<observables;i++){
    for(j=0;j<repetitions;j++){
      mean(i) += matrix(j,i);
    }
    mean(i) /= (double) repetitions;
  }
  for(i=0;i<observables;i++){
    xi_mean = mean(i);
    for(j=0;j<observables;j++){
      xj_mean = mean(j);
      for(k=0;k<repetitions;k++){
	covariance(i,j) += (matrix(k,i) - xi_mean)*(matrix(k,j) - xj_mean);
      }
      covariance(i,j) /= (double) repetitions;
    }
  }
  for(i=0;i<observables;i++){
    xi_mean = mean(i);
    xi_std = sqrt(covariance(i,i));
    for(j=0;j<repetitions;j++){
      tilde(j,i) = (matrix(j,i) - xi_mean)/xi_std;
    }
  }
}
void covariance_function(const arma::mat &matrix, arma::mat &covariance){
  int repetitions = matrix.n_rows;
  int observables = matrix.n_cols;
  covariance = arma::zeros<arma::mat>(observables,observables);

  for(int i=0;i<observables;i++){
    for(int j=0;j<observables;j++){
      for(int k=0;k<repetitions;k++){
	covariance(i,j) += matrix(k,i)*matrix(k,j);
      }
      covariance(i,j) /= (double) repetitions;
    }
  }
}
void sort_eigen_function(arma::vec &eigval, arma::mat &eigvec){
  int i, j, sort,
    observables = eigval.n_elem;
  arma::mat eigsort = arma::zeros<arma::mat>(observables,observables);
  arma::mat eigval_matrix = arma::zeros<arma::mat>(observables,2);

  for(i=0;i<observables;i++){
    eigval_matrix(i,0) = eigval(i);
    eigval_matrix(i,1) = i;
  }
  eigval_matrix = arma::sort(eigval_matrix, "descending", 0);
  for(i=0;i<observables;i++){
    sort = eigval_matrix(i,1);
    for(j=0;j<observables;j++){
      eigsort(j,i) = eigvec(j,sort);
    }
    eigval(i) = eigval_matrix(i,0);
  }
  for(i=0;i<observables;i++){
    for(j=0;j<observables;j++){
      eigvec(i,j) = eigsort(i,j);
    }
  }
}
// END: Functions to create the y~ matrix and conduct PCA


/******

       Functions for statistics

 ******/

double median_width_function(const arma::mat &function, double factor){
  int points = function.n_rows - 1;
  double median_width,sum=0.0;
  arma::vec tally = arma::zeros<arma::vec>(points);
  for(int i=0;i<points;i++){
    tally(i) = (function(i+1,1) + function(i,1))*(function(i+1,0) - function(i,0))/2.0;
    sum += tally(i);
    tally(i) = sum;
  }
  sum *= factor;
  if(tally(points-1) <= 0.0){
    printf("median_Width_function has been given a function with area <= 0.0. Terminating...\n");
    exit(1);
  }
  else{
    for(int i=0;i<points;i++){
      if(tally(i)>sum){
	median_width = (function(i+1,0) + function(i,0))/2.0;
	break;
      }
    }
  }
  return median_width;
}
double median_width_fabs_function(const arma::mat &function, double factor){
  int points = function.n_rows - 1;
  double median_width,sum=0.0;
  arma::vec tally = arma::zeros<arma::vec>(points);
  for(int i=0;i<points;i++){
    tally(i) = fabs((function(i+1,1) + function(i,1)))*(function(i+1,0) - function(i,0))/2.0;
    sum += tally(i);
    tally(i) = sum;
  }
  sum *= factor;
  if(tally(points-1) <= 0.0){
    printf("median_Width_function has been given a function with area <= 0.0. Terminating...\n");
    exit(1);
  }
  else{
    for(int i=0;i<points;i++){
      if(tally(i)>sum){
	median_width = (function(i+1,0) + function(i,0))/2.0;
	break;
      }
    }
  }
  return median_width;
}
double zeroth_moment(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, f, dx;
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  return zero;
}
double first_moment(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, first = 0.0, f, x, dx;
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    first += f*x*dx;
  }
  first /= zero;
  return first;
}
double second_moment(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, first = 0.0, second = 0.0;
  double f, x, dx;
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    first += f*x*dx;
  }
  first /= zero;
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    second += (x - first)*(x - first)*f*dx;
  }
  second /= zero;
  return second;
}
// END: Functions for statistics

