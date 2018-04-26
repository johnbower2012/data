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
arma::mat kernal_sqare_exp_function(arma::mat A, arma::mat B, arma::vec hyperp_vec);
arma::mat kernal_derivative_function(arma::mat kernal_func(arma::mat, arma::mat, arma::vec), arma::mat A, arma::mat B, arma::vec hyperp_vec, int selection, double precision);

/**GAUSSIAN PROCESS FUNCTIONS
**********************************************/
arma::mat gaussian_process_solver_regression(
	arma::mat kernal_func(arma::mat, arma::mat, arma::vec), arma::mat X_mat, arma::mat Y_mat, arma::mat X_s_mat, 
	int samples, arma::vec hyperp_vec, double epsilon, 
	arma::mat regression_func(arma::mat, arma::vec), int power);
arma::mat gaussian_process_solver_basic(arma::mat X_mat, arma::mat Y_mat, arma::mat X_s_mat, int samples, arma::vec hyperp_vec, double epsilon);

/**REGRESSION FUNCTIONS
**********************************************/
arma::mat regression_polynomial_function(arma::mat X_mat, arma::vec parameters);

/**WRITE TO FILE FUNCTIONS
**********************************************/
void write_output(arma::mat output_mat, int test, int param, int observables, int samples, std::string outfilename);

/**SEARCH FUNCTIONS
**********************************************/
int find_max(arma::vec input_vec);
int find_zero(arma::vec input_vec);

/**LIKELIHOOD FUNCTIONS
**********************************************/
arma::vec log_likelihood_function(arma::mat Y_mat, arma::mat X_mat, int resolution, arma::vec hyperp_vec, arma::mat hyperp_range_vec, int selection);
arma::vec log_likelihood_derivative(arma::mat Y_mat, arma::mat X_mat, int resolution, arma::vec hyperp_vec, arma::mat hyperp_range_vec, int selection, double precision);

#endif
