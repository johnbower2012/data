#ifndef GAUSSIANPROCESS_H
#define GAUSSIANPROCESS_H

#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include "time.h"
#include "armadillo"

#include "gaussianprocess.hpp"


/**KERNAL FUCNTIONS
**********************************************/
arma::mat kernel_sqare_exp_function(arma::mat A, arma::mat B, arma::vec hyperp_vec);
arma::mat kernel_sqare_exp_noise_function(arma::mat A, arma::mat B, arma::vec hyperp_vec);
arma::mat kernel_derivative_function(arma::mat kernal_func(arma::mat, arma::mat, arma::vec), arma::mat A, arma::mat B, arma::vec hyperp_vec, int selection, double precision);

/**GAUSSIAN PROCESS FUNCTIONS
**********************************************/
arma::mat gaussian_process_solver_regression(
	arma::mat kernal_func(arma::mat, arma::mat, arma::vec), arma::mat X_mat, arma::mat Y_mat, arma::mat X_s_mat, 
	int samples, arma::vec hyperp_vec, double epsilon, 
	arma::mat regression_func(arma::mat, arma::mat), arma::mat beta);
arma::mat gaussian_process_solver_basic(arma::mat X_mat, arma::mat Y_mat, arma::mat X_s_mat, int samples, arma::vec hyperp_vec, double epsilon);

/**REGRESSION FUNCTIONS
**********************************************/
arma::mat regression_polynomial_function(arma::mat X_mat, arma::vec parameters);
arma::mat regression_linear_function(arma::mat X_mat, arma::mat beta);

/**WRITE TO FILE FUNCTIONS
**********************************************/
void write_output(arma::mat output_mat, int test, int param, int observables, int samples, std::string outfilename);
void write_trainset(arma::mat xvec_train_mat, arma::mat yvec_train_mat, std::string outfilename);

/**SEARCH FUNCTIONS
**********************************************/
int find_max(arma::vec input_vec);
arma::vec find_zeros(arma::vec input_vec);

#endif
