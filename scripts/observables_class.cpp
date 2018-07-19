#include<iostream>
#include<fstream>
#include<cmath>
#include "armadillo"

void tilde_function(const arma::mat &matrix, const arma::vec &error, arma::vec &mean, arma::mat &tilde);
void tilde_function(const arma::mat &matrix, arma::mat &covariance, arma::vec &mean, arma::mat &tilde);
void covariance_function(const arma::mat &matrix, arma::mat &covariance);
void sort_eigen_function(arma::vec &eigval, arma::mat &eigvec);

double zeroth_moment(const arma::mat &function);
double first_moment(const arma::mat &function);
double second_moment(const arma::mat &function);
void obs_matrix_moments(int files, int obs_file, arma::mat *&val_matrix, const arma::vec &delY_vec, arma::mat &obs_matrix);

double zeroth_moment_fabs(const arma::mat &function);
double first_moment_fabs(const arma::mat &function);
double second_moment_fabs(const arma::mat &function);
void obs_matrix_moments_fabs(int files, int obs_file, arma::mat *&val_matrix, const arma::vec &delY_vec, arma::mat &obs_matrix);

arma::vec average_columns(arma::mat input);
void load_file(int files, int lines, int runs, std::string *&infilename, arma::vec &delY_vec, arma::mat *&val_matrix);
void print_file(std::string outfilename, std::string title, arma::vec vector);
void print_file(std::string outfilename, arma::vec vector);
void print_file(std::string outfilename, std::string title, arma::mat matrix);
void print_file(std::string outfilename, arma::mat matrix);
void print_fractional_sum(arma::vec vector);


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
      tilde(i,j) = (matrix(i,j) - mean(j))/error(j);
    }
  }
}
void tilde_function_input(const arma::mat &matrix, const arma::vec &error, arma::vec &mean, arma::mat &tilde){
  int repetitions = matrix.n_rows;
  int observables = matrix.n_cols;
  tilde = arma::zeros<arma::mat>(repetitions,observables);

  for(int i=0;i<repetitions;i++){
    for(int j=0;j<observables;j++){
      tilde(i,j) = (matrix(i,j) - mean(j))/error(j);
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

double zeroth_moment(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, f, dx;
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  return zero;
}
double first_moment(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, first = 0.0, f, x, dx;

  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }

  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
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
    f = (function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }

  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    first += f*x*dx;
  }
  first /= zero;
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    second += (x - first)*(x - first)*f*dx;
  }
  second /= zero;
  return second;
}
void obs_matrix_moments(int files, int obs_file, arma::mat *&val_matrix, const arma::vec &delY_vec, arma::mat &obs_matrix){
  int i,j,k;
  int runs = val_matrix[0].n_rows;
  int lines = delY_vec.n_elem;
  arma::mat function = arma::zeros<arma::mat>(lines,2);
  obs_matrix = arma::zeros<arma::mat>(runs,obs_file*files);
  for(i=0;i<files;i++){
    for(j=0;j<runs;j++){
      for(k=0;k<lines;k++){
	function(k,0) = delY_vec(k);
	function(k,1) = val_matrix[i](j,k);
      }
      obs_matrix(j,i*obs_file) = zeroth_moment(function);
      obs_matrix(j,i*obs_file+1) = first_moment(function);
      obs_matrix(j,i*obs_file+2) = second_moment(function);
    }
  }  
}
//ABS VERSIONS
double zeroth_moment_fabs(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, f, dx;
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  return zero;
}
double first_moment_fabs(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, first = 0.0, f, x, dx;

  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }

  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    first += f*x*dx;
  }
  //  first /= zero;
  return first;
}
double second_moment_fabs(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, first = 0.0, second = 0.0;
  double f, x, dx;
  
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    first += f*x*dx;
  }
  //  first /= zero;
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    second += (x - first)*(x - first)*f*dx;
  }
  //  second /= zero;
  return second;
}
void obs_matrix_moments_fabs(int files, int obs_file, arma::mat *&val_matrix, const arma::vec &delY_vec, arma::mat &obs_matrix){
  int i,j,k;
  int runs = val_matrix[0].n_rows;
  int lines = delY_vec.n_elem;
  arma::mat function = arma::zeros<arma::mat>(lines,2);
  obs_matrix = arma::zeros<arma::mat>(runs,obs_file*files);
  for(i=0;i<files;i++){
    for(j=0;j<runs;j++){
      for(k=0;k<lines;k++){
	function(k,0) = delY_vec(k);
	function(k,1) = val_matrix[i](j,k);
      }
      obs_matrix(j,i*obs_file) = zeroth_moment_fabs(function);
      obs_matrix(j,i*obs_file+1) = first_moment_fabs(function);
      obs_matrix(j,i*obs_file+2) = second_moment_fabs(function);
    }
  }  
}

// END: Functions for statistics
arma::vec average_columns(arma::mat input){
  int rows = input.n_rows;
  int columns = input.n_cols;
  arma::vec avg = arma::zeros<arma::vec>(columns);
  for(int i=0;i<columns;i++){
    for(int j=0;j<rows;j++){
      avg(i) += input(j,i);
    }
    avg(i) /= (double) rows;
  }
  return avg;
}
void load_file(int files, int lines, int runs, std::string *&infilename, arma::vec &delY_vec, arma::mat *&val_matrix){
  std::ifstream ifile;
  if(val_matrix!=nullptr){
    delete[] val_matrix;
  }
  val_matrix = new arma::mat[files];
  int i, j, k;

  for(i=0;i<files;i++){
    val_matrix[i] = arma::zeros<arma::mat>(runs,lines);
    ifile.open(infilename[i]);
    printf("+++++ %s LOADED +++++\n", infilename[i].c_str());
    for(j=0;j<lines;j++){
      for(k=0;k<runs+1;k++){
	if(k==0){
	  ifile >> delY_vec(j);
	} else{
	  ifile >> val_matrix[i](k-1,j);
	}
      }
    }
    ifile.close();
  }
}
void print_file(std::string outfilename, std::string title, arma::vec vector){
  int elem = vector.n_elem;
  std::ofstream ofile;
  ofile.open(outfilename);
  ofile << title.c_str() << '\n';
  for(int i=0;i<elem;i++){
    ofile << ' ' << vector(i) << '\n';
  }
  ofile.close();
}
void print_file(std::string outfilename, arma::vec vector){
  int elem = vector.n_elem;
  std::ofstream ofile;
  ofile.open(outfilename);
  for(int i=0;i<elem;i++){
    ofile << ' ' << vector(i) << '\n';
  }
  ofile.close();
}
void print_file(std::string outfilename, std::string title, arma::mat matrix){
  int rows = matrix.n_rows;
  int cols = matrix.n_cols;
  std::ofstream ofile;
  ofile.open(outfilename);
  ofile << title.c_str() << '\n';
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      ofile << ' ' << matrix(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
}
void print_file(std::string outfilename, arma::mat matrix){
  int rows = matrix.n_rows;
  int cols = matrix.n_cols;
  std::ofstream ofile;
  ofile.open(outfilename);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      ofile << ' ' << matrix(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
}

//print running fractional total of vector
void print_fractional_sum(arma::vec vector){
  int length = vector.n_elem;
  arma::vec sum = arma::zeros<arma::vec>(length);
  for(int i=0;i<length;i++){
    sum(i) = vector(i);
    sum(i) += sum(abs(i-1));
  }
  for(int i=0;i<length;i++){
    sum(i) /= sum(length-1);
   }
  sum.print("+++++ eigval fractional sum +++++");
}
