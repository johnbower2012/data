#ifndef PCA_H
#define PCA_H

#include<iomanip>
#include<cstdlib>
#include<random>
#include<chrono>
#include<fstream>
#include<string>
#include "armadillo"

//creates vector of mean values for each column of a matrix
arma::vec calculate_mean_function(arma::mat val_matrix);
//creates covariance matrix for a matrix of n-repitions of p-dim row vectors
arma::mat calculate_covariance_function(arma::mat val_matrix, arma::vec mean_vec);
//calculate zmatrix, with mean zero and unit std, i.e. [z = (y-ybar)/ystd]
arma::mat calculate_zmatrix_function(arma::mat val_matrix, arma::vec mean_vec, arma::mat cov_matrix);
//Calculate valmatrix from zmatrix
arma::mat calculate_valmatrix_function(arma::mat zval_matrix, arma::vec mean_vec, arma::mat cov_matrix);
//creates covariance matrix for a matrix of n-repitions of p-dim row vectors, with zero column-wise mean
arma::mat calculate_zcovariance_function(arma::mat zval_matrix);
//sort eigenvectors according to eigenvalues, biggest to smallest
arma::mat sort_eigenvectors_function(arma::vec eigval_vec, arma::mat eigvec_matrix);
//full transformation from val_matrix to tval_matrix
arma::mat calculate_tmatrix_function(arma::mat val_matrix);
//Transformation from zval_matrix to tval_matrix
arma::mat calculate_ztotmatrix_function(arma::mat zval_matrix);
//full transformation from tval_matrix to val_matrix
arma::mat calculate_ttovalmatrix_function(arma::mat tval_matrix, arma::mat eigvec_matrix, arma::vec mean_vec, arma::mat cov_matrix);



#endif
