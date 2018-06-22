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
double zeroth_moment_fabs(const arma::mat &function);
double first_moment_fabs(const arma::mat &function);
double second_moment_fabs(const arma::mat &function);

void obs_matrix_median_widths(int files, int obs_file, arma::mat *&val_matrix, const arma::vec &delY_vec, arma::mat &obs_matrix);
void obs_matrix_moments(int files, int obs_file, arma::mat *&val_matrix, const arma::vec &delY_vec, arma::mat &obs_matrix);
void obs_matrix_moments_fabs(int files, int obs_file, arma::mat *&val_matrix, const arma::vec &delY_vec, arma::mat &obs_matrix);

arma::vec average_columns(arma::mat input);
void load_file(int files, int lines, int runs, std::string *&infilename, arma::vec &delY_vec, arma::mat *&val_matrix);
void print_file(std::string outfilename, std::string title, arma::mat matrix);
void print_file(std::string outfilename, arma::mat matrix);
void print_fractional_sum(arma::vec vector);

int main(int argc, char* argv[]){
  std::string *infilename,*outfilename;
  int observables,obs_file,files,lines,runs;
  int i,j,k;
  infilename = new std::string[observables];
  outfilename = new std::string[observables];
  std::string dest_folder;

  files=4;
  obs_file=3;
  observables=files*obs_file;

  if(argc<4+files){
    std::cout << "Improper input. Enter also 'lines runs observables ifn*[7]' on same line." << std::endl;
    exit(1);
  }
  else{
    infilename = new std::string[files];
    dest_folder = argv[1];
    lines=atoi(argv[2]);
    runs=atoi(argv[3]);
    files=atoi(argv[4]);
    for(i=0;i<files;i++){
      infilename[i]=argv[5+i];
    }
  }
  printf("Arguments read in as:\n");
  printf("observables.x %s %d %d %d",dest_folder.c_str(),lines,runs,files);
  for(i=0;i<files;i++){
    printf(" %s",infilename[i].c_str());
  }
  printf("\n");

  /*********
	     LOAD FILE
	     CALCULATE OBS_MAT
  *********/
  arma::mat *val_matrix, obs_matrix;
  arma::vec delY_vec, obs_error;
  val_matrix = new arma::mat[files];
  obs_matrix = arma::zeros<arma::mat>(runs,observables);
  delY_vec = arma::zeros<arma::vec>(lines);
  obs_error = arma::zeros<arma::vec>(observables);
  //all observables that aren't pipi, reduce unc by factor of ten.

  //pipi
  int part=0;
  obs_error(part*obs_file) = 0.01;
  obs_error(1+part*obs_file) = 0.01;
  obs_error(2+part*obs_file) = 0.005;
  //ppbar
  part=1;
  obs_error(part*obs_file) = 0.005;
  obs_error(1+part*obs_file) = 0.005;
  obs_error(2+part*obs_file) = 0.001;
  //pK
  part=2;
  obs_error(part*obs_file) = 0.001;
  obs_error(1+part*obs_file) = 0.001;
  obs_error(2+part*obs_file) = 0.0005;
  //KK
  part=3;
  obs_error(part*obs_file) = 0.005;
  obs_error(1+part*obs_file) = 0.005;
  obs_error(2+part*obs_file) = 0.001;


  load_file(files, lines, runs, infilename, delY_vec, val_matrix);
  obs_matrix_moments_fabs(files,obs_file,val_matrix,delY_vec,obs_matrix);

  arma::vec avg_mom = average_columns(obs_matrix);
  avg_mom.print("avg_mom");
  obs_error = 0.1*avg_mom;
  obs_error.print("obs");



  /*********
	     CONDUCT PCA
	     PRINT
  *********/

 
  arma::mat tval_matrix,zcov_matrix,eigvec_matrix,print_matrix;
  arma::vec eigval_vec,mean_vec;
  print_matrix = arma::zeros<arma::mat>(observables,observables+1);
  tval_matrix = arma::zeros<arma::mat>(observables,observables);
  zcov_matrix = arma::zeros<arma::mat>(observables,observables);
  eigvec_matrix = arma::zeros<arma::mat>(observables,observables);
  eigval_vec = arma::zeros<arma::vec>(observables);
  mean_vec = arma::zeros<arma::vec>(observables);

  tval_matrix = calculate_tmatrix_function(obs_matrix, obs_error, eigval_vec, eigvec_matrix, mean_vec, zcov_matrix);

  print_fractional_sum(eigval_vec);
  eigval_vec.print("+++++ eigvalues +++++");
  eigvec_matrix.print("+++++ eigvectors +++++");

  for(i=0;i<observables;i++){
    print_matrix(i,0) = eigval_vec(i);
    for(j=0;j<observables;j++){
      print_matrix(i,j+1) = eigvec_matrix(j,i);
    }
  }
  std::string printname;
  std::string title;
  printname = "moments_model_data.dat";
  title = "#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,title,obs_matrix);

  printname = "moments_model_eigvec.dat";
  title = "#colvec\n#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,title,eigvec_matrix);

  printname = "moments_model_pca.dat";
  title = "#eigval col -- eigvec rows\n#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,title,print_matrix);

  print_matrix = obs_matrix*eigvec_matrix;
  printname = "moments_model_z.dat";
  title = "#moments_matrix*eigenvectors_matrix";
  print_file(printname,title,print_matrix);

  delete[] val_matrix;
  delete[] infilename;
  delete[] outfilename;

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
    f = (function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  return zero;
}
double first_moment(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, first = 0.0, f, x, dx;
  /*
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  */
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    first += f*x*dx;
  }
  //  first /= zero;
  return first;
}
double second_moment(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, first = 0.0, second = 0.0;
  double f, x, dx;
  /*
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  */
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    first += f*x*dx;
  }
  //first /= zero;
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    second += (x - first)*(x - first)*f*dx;
  }
  //second /= zero;
  return second;
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
  /*
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  */
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    first += f*x*dx;
  }
  //first /= zero;
  return first;
}
double second_moment_fabs(const arma::mat &function){
  int points = function.n_rows - 1;
  double zero = 0.0, first = 0.0, second = 0.0;
  double f, x, dx;
  /*
  for(int i=0;i<points;i++){
    f = fabs(function(i+1,1) + function(i,1))/2.0;
    dx = function(i+1,0) - function(i,0);
    zero += f*dx;
  }
  */
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    first += f*x*dx;
  }
  //first /= zero;
  for(int i=0;i<points;i++){
    f = (function(i+1,1) + function(i,1))/2.0;
    x = (function(i+1,0) + function(i,0))/2.0;
    dx = function(i+1,0) - function(i,0);
    second += (x - first)*(x - first)*f*dx;
  }
  //second /= zero;
  return second;
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

void obs_matrix_median_widths(int files, int obs_file, arma::mat *&val_matrix, const arma::vec &delY_vec, arma::mat &obs_matrix){
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
      obs_matrix(j,i*obs_file) = median_width_fabs_function(function,0.5);
      obs_matrix(j,i*obs_file+1) =  median_width_fabs_function(function,2.0/3.0) - median_width_fabs_function(function,1.0/3.0);
    }
  }
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
