/*
 *
 *
	This file uses "pca.cpp" and "armadillo" to conduct Principle Component Analysis. PCA is
	a rotation whereby the new axes are in the directions of greatest variance for the data,
	ordered {t_i} such that t_1 is in the direction of greatest variance, t_2 the next, and
	so forth.

	One begins with a matrix of dimensions n-by-d, here called the "val_matrix," where there 
	are n-repitions of a d-dimensional row vector, that is the 'i'th row is a different 
	measurement of a set of d values, {x_j}_i, and the 'j'th column is a set of n values, 
	{x_i}_j, each x_i a new measurement of the 'j'th dimension.

	PCA functions by first translating and scaling each column to have zero mean and unit 
	variance, here called the "zval_matrix". Then, it is desired to find weights, W = {w_i},
	each |w_i| = 1, such that T = Z*W has as its axes the directions of great variance. In 
	this aim, the covariance matrix for the zval_matrix is constructed and its eigenvalues 
	& vectos found. If w_1 is the eigenvector corresponding to lambda_1 = max{lambda_i}, 
	then t_1 will be the PC in the direction of greatest variance. If we also then arrange 
	lambda_1 > lambda_2 > ... > lambda_n, with corresponding {w_i}, we then have our {t_i} 
	vectors in order of greatest to least variance.

	The mathematical proof of this arrangement follows from the fact that we want our t_1
	to satisfy t_1 = max{|w_T*Z_T*Z*w|} with |w| = 1, which is the equivalent of 
	t_1 = max{|w_T*Z_T*Z*w/(w_T*w)|} since |w| = 1. This can be recognized as the Rayleigh
	Quotient, where the maximum is known to be the largest eigenvalue of Z_T*Z and occurs
	when w is the corresponding eigenvector. Lastly, Z_T*Z = SUM_over_k(z_ki*z_kj) = cov(Z),
	and so the above has been shown (mostly).
 *
 *
 */

#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<random>
#include<chrono>
#include<fstream>
#include<string>
#include "armadillo"
#include "pricomana.cpp"

int main(int argc, char* argv[]){
  int reps, dims,observables,lines,runs;
  int i,j,k;
  std::string *infilename, outfilename;
  std::ofstream ofile;
  std::ifstream ifile;

  observables=7;
  
  if(argc<5+observables){
    std::cout << "Improper input. Enter also 'lines runs observables ifn*[7] ofn' on same line."  << std::endl;
    exit(1);
  }
  else{
    infilename = new std::string[observables];
    lines=atoi(argv[1]);
    runs=atoi(argv[2]);
    observables=atoi(argv[3]);
    for(i=0;i<observables;i++){
      infilename[i] = argv[4+i];
    }
    outfilename = argv[4+observables];
  }
  //TEST LESS RUNS AND LINES
  lines=21;
  //runs=999;
  observables=1;

  //initialize required vectors and matrices
  arma::mat *val_matrix,*tval_matrix,*cov_matrix,*eigvec_matrix;
  arma::vec *eigval_vec,*mean_vec;
  arma::mat print_matrix = arma::zeros<arma::mat>(lines,lines+1);
  arma::vec delY = arma::zeros<arma::vec>(lines);
  val_matrix = new arma::mat[observables];
  tval_matrix = new arma::mat[observables];
  cov_matrix = new arma::mat[observables];
  eigvec_matrix = new arma::mat[observables];
  eigval_vec = new arma::vec[observables];
  mean_vec = new arma::vec[observables];
  
  for(i=0;i<observables;i++){
    val_matrix[i] = arma::zeros<arma::mat>(runs,lines);
    tval_matrix[i] = arma::zeros<arma::mat>(runs,lines);
    cov_matrix[i] = arma::zeros<arma::mat>(lines,lines);
    eigvec_matrix[i] = arma::zeros<arma::mat>(lines,lines);
    eigval_vec[i] = arma::zeros<arma::vec>(lines);
    mean_vec[i] = arma::zeros<arma::vec>(lines);
  }

  //load data from file
  for(i=0;i<observables;i++){
    ifile.open(infilename[i]);
    for(j=0;j<lines;j++){
      for(k=0;k<runs+1;k++){
	if(k==0){
	  ifile >> delY(j);
	}
	else{
	  ifile >> val_matrix[i](k-1,j);
	}
      }
    }
    ifile.close();
  }
  
  //calculate the PCA Transformed matrix
  for(i=0;i<observables;i++){
    printf("+++++ %s +++++\n",infilename[i].c_str());
    tval_matrix[i] = calculate_tmatrix_function(val_matrix[i], eigval_vec[i], eigvec_matrix[i], mean_vec[i], cov_matrix[i]);
    eigval_vec[i].print();
  }

  //write output to file
  std::string delimiter;
  std::size_t position;
  
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

  delete[] val_matrix;
  delete[] tval_matrix;
  delete[] eigval_vec;
  delete[] eigvec_matrix;
  delete[] mean_vec;
  delete[] cov_matrix;
  delete[] infilename;
  
  return 0;

}
