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
arma::mat kernel_periodic_function(arma::mat A, arma::mat B, arma::vec hyperp_vec);
arma::mat kernel_periodic_decayaway_function(arma::mat A, arma::mat B, arma::vec hyperp_vec);
arma::mat kernel_rational_quadratic_function(arma::mat A, arma::mat B, arma::vec hyperp_vec);
arma::mat kernel_co2_function(arma::mat A, arma::mat B, arma::vec hyperp_vec);
arma::mat kernel_derivative_function(arma::mat kernal_func(arma::mat, arma::mat, arma::vec), arma::mat A, arma::mat B, arma::vec hyperp_vec, int selection, double precision);

/**GAUSSIAN PROCESS FUNCTIONS
**********************************************/
arma::mat gaussian_process_solver_regression(
	arma::mat kernal_func(arma::mat, arma::mat, arma::vec), arma::mat X_mat, arma::mat Y_mat, arma::mat X_s_mat, 
	int samples, arma::vec hyperp_vec, double epsilon, 
	arma::mat regression_func(arma::mat, arma::vec), int power);
arma::mat gaussian_process_solver_basic(arma::mat X_mat, arma::mat Y_mat, arma::mat X_s_mat, int samples, arma::vec hyperp_vec, double epsilon);
void gaussian_process(const arma::mat &x_mat, const arma::vec &y_vec, arma::mat kernel_function(arma::mat,arma::mat,arma::vec), arma::vec hyperp_param, double sigma_n, const arma::mat &x_star_mat,
		      arma::vec &mean, arma::mat &variance, double &log_likelihood);

/**REGRESSION FUNCTIONS
**********************************************/
arma::mat regression_polynomial_function(arma::mat X_mat, arma::vec parameters);

/**WRITE TO FILE FUNCTIONS
**********************************************/
void write_output(arma::mat output_mat, int test, int param, int observables, int samples, std::string outfilename);
void write_trainset(arma::mat xvec_train_mat, arma::mat yvec_train_mat, std::string outfilename);

/**SEARCH FUNCTIONS
**********************************************/
int find_max(arma::vec input_vec);
arma::vec find_zeros(arma::vec input_vec);

/**LIKELIHOOD FUNCTIONS
**********************************************/
long double log_likelihood_function(arma::mat Y_mat, arma::mat X_mat, arma::mat kernal_func(arma::mat, arma::mat, arma::vec), arma::vec hyperp_vec);
arma::vec log_likelihood_function(arma::mat Y_mat, arma::mat X_mat, int resolution, arma::mat kernal_func(arma::mat, arma::mat, arma::vec), arma::vec hyperp_vec, arma::mat hyperp_range_vec, int selection);
arma::vec log_likelihood_spec_function(arma::mat Y_mat, arma::mat X_mat, arma::mat kernal_func(arma::mat, arma::mat, arma::vec), arma::vec hyperp_vec, arma::vec hyperp_position_vec, int selection);
arma::vec log_likelihood_derivative(arma::mat Y_mat, arma::mat X_mat, int resolution, arma::vec hyperp_vec, arma::mat hyperp_range_vec, int selection, double precision);

#endif
