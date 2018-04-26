#include<iostream>
#include<fstream>
#include<cmath>
#include "armadillo"
#include "pricomana.cpp"

void tilde_function(const arma::mat &matrix, const arma::vec &error, arma::vec &mean, arma::mat &tilde);
void tilde_function(const arma::mat &matrix, arma::mat &covariance, arma::vec &mean, arma::mat &tilde);
void covariance_function(const arma::mat &matrix, arma::mat &covariance);
void sort_eigen_function(arma::vec &eigval, arma::mat &eigvec);

double median_width_function(const arma::mat &function, double factor);
double median_width_fabs_function(const arma::mat &function, double factor);
double zeroth_moment(const arma::mat &function);
double first_moment(const arma::mat &function);
double second_moment(const arma::mat &function);

void load_file(int files, int lines, int runs, std::string *&infilename, arma::vec &delY_vec, arma::mat *&val_matrix);

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
  arma::mat *val_matrix, *obs_matrix;
  arma::vec delY_vec;
  val_matrix = new arma::mat[observables];
  obs_matrix = new arma::mat[observables];
  delY_vec = arma::zeros<arma::vec>(lines);
  load_file(observables, lines, runs, infilename, delY_vec, val_matrix);

  for(i=0;i<observables;i++){
    for(j=0;j<5;j++){
      for(k=0;k<10;k++){
	printf("%f ",val_matrix[i](j,k));
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("Pass\n");
  delY_vec.print();
  /*
    for(i=0;i<observables;i++){
    val_matrix[i] = arma::zeros<arma::mat>(runs,lines);
    obs_matrix[i] = arma::zeros<arma::mat>(runs,2);
    ifile.open(infilename[i]);
    printf("+++++ %s LOADED +++++\n", infilename[i].c_str());
    for(j=0;j<lines;j++){
      for(k=0;k<runs;k++){
	if(k==0){
	  ifile >> delY_vec(j);
	} else{
	  ifile >> val_matrix[i](k-1,j);
	}
      }
    }
    ifile.close();
  }
  */
  
  arma::mat function = arma::zeros<arma::mat>(runs,2);
  for(i=0;i<observables;i++){
    for(j=0;j<runs;j++){
      for(k=0;k<lines;k++){
	function(k,0) = delY_vec(k);
	function(k,1) = val_matrix[i](j,k);
      }
      obs_matrix[i](j,0) = median_width_fabs_function(function,0.5);
      obs_matrix[i](j,1) =  median_width_fabs_function(function,2.0/3.0) - median_width_fabs_function(function,1.0/3.0);
    }
  }

  arma::mat *tval_matrix,*cov_matrix,*eigvec_matrix;
  arma::vec *eigval_vec,*mean_vec;
  arma::mat print_matrix = arma::zeros<arma::mat>(2,2);
  tval_matrix = new arma::mat[observables];
  cov_matrix = new arma::mat[observables];
  eigvec_matrix = new arma::mat[observables];
  eigval_vec = new arma::vec[observables];
  mean_vec = new arma::vec[observables];
  for(i=0;i<observables;i++){
    tval_matrix[i] = arma::zeros<arma::mat>(2,2);
    cov_matrix[i] = arma::zeros<arma::mat>(2,2);
    eigvec_matrix[i] = arma::zeros<arma::mat>(2,2);
    eigval_vec[i] = arma::zeros<arma::vec>(2);
    mean_vec[i] = arma::zeros<arma::vec>(2);
  }
  for(i=0;i<observables;i++){
    printf("+++++ %s eigvalues +++++\n",infilename[i].c_str());
    tval_matrix[i] = calculate_tmatrix_function(obs_matrix[i], eigval_vec[i], eigvec_matrix[i], mean_vec[i], cov_matrix[i]);
    eigval_vec[i].print();
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
  delete[] obs_matrix;
  delete[] infilename;
  delete[] outfilename;
  delete[] tval_matrix;
  delete[] eigval_vec;
  delete[] eigvec_matrix;
  delete[] mean_vec;
  delete[] cov_matrix;
  
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
