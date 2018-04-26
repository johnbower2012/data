/* I just realized there is a fairly large flaw in my output
	matrix in the GPS. I account for multiple obervables in 
	the mean and std output, but not the generated functions.
	Update in future. Easy fix.

	***********************************************************/


#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include<fstream>
#include<string>
#include "time.h"
#include "armadillo"
#include "gaussianprocess.cpp"

std::ofstream ofile;

/*
double function(arma::vec param){
	double value = param(1)*param(1)*sin(1.5*param(0));
	return value;
}
*/
double function(arma::vec param){
	double value = param(1) + sin(1.5*param(0));
	return value;
}
/*
double function(arma::vec param){
	double value = param(0) + sin(1.5*param(0));
	return value;
}
*/



int main(int argc, char* argv[]){
	//Variables for use throughout the program
	int i, j, k,
		train, test, 
		param, observables, hyperp, power, basic_lr,
		samples;
	unsigned seed;
	double 	x, x_ini, x_fin, dx,
			sigma_f, sigma_n, l, epsilon;

	//For writing to file
	std::string outfilename;

	//For timing, if desired
	clock_t start, finish;

	//Test to ensure the proper input. Terminate program if command line inputs are not given
	if(argc<5){
		std::cout << "Improper entry. Please also enter 'train test x_ini x_fin' on same line." << std::endl;
		exit(1);
	}
	else{
		train = atoi(argv[1]);
		if(train<2){
			std::cout << "Improper entry. 'train' must be greater than 1." << std::endl;
			exit(1);
		}
		else{
			//manual user input
			train = atoi(argv[1]);
			test = atoi(argv[2]);
			x_ini = atof(argv[3]);
			x_fin = atof(argv[4]);
			//power = atoi(argv[5]);
			power = 0;
		}
	}

	/*SYSTEM and OUTPUT--
		specify--
			input dimensionality: 				param
			output dimensionality: 				observables
			sets of generated output values:	samples

			numerical stability addition:		epsilon

		generated from computer clock--
			ranomd seed for N(0,1): 			seed
	***********************************************************/

		samples = 3;
		param = 2;
		observables = 1;

		epsilon = 1e-8;

		seed = std::chrono::system_clock::now().time_since_epoch().count();

	/*HYPERPARAMETERS:
		specify--
			number of hyperparameters:			hyperp
			noise amplitude:					sigma_n
			variation amplitude:				sigma_f
			length scale of variation:			l;			
	***********************************************************/

		hyperp = 3;
		sigma_n = 0.001;
		sigma_f = 1;
		l = 1.0;

	//random number generator
	std::default_random_engine generator(seed);
	std::normal_distribution<double> dist(0,1.0);
	std::uniform_real_distribution<double> dist_r(x_ini,x_fin);

	//armadillo matrices and vectors for computational use
	arma::mat 	xvec_train_mat = arma::zeros<arma::mat>(train,param),
				yvec_train_mat = arma::zeros<arma::mat>(train,observables),
				xvec_test_mat = arma::zeros<arma::mat>(test,param),
				yvec_test_mat = arma::zeros<arma::mat>(test,observables),

				output_mat = arma::zeros<arma::mat>(test,param + 3 + samples);

	arma::vec	hyperp_vec = arma::zeros<arma::vec> (hyperp),
				param_x_vec = arma::zeros<arma::vec> (param);

	//hyperparamter vector for use in 'function'
	hyperp_vec(0) = sigma_f;
	hyperp_vec(1) = l;
	hyperp_vec(2) = sigma_n;

	//Generate Functions
	dx = x_fin - x_ini;
	for(i=0;i<train;i++){
		for(j=0;j<param;j++){
			xvec_train_mat(i,j) = /*x_ini + dx*(double) i/(double) (train-1);//*/dist_r(generator);
		}
	}
	for(i=0;i<test;i++){
		for(j=0;j<param;j++){
			xvec_test_mat(i,j) = /*x_ini + dx*(double) i/(double) (test-1);//*/dist_r(generator);
		}
	}
	for(i=0;i<train;i++){
		for(j=0;j<param;j++){
			param_x_vec(j) = xvec_train_mat(i,j);
		}
		yvec_train_mat(i,0) = function(param_x_vec);
	}
	for(i=0;i<test;i++){
		for(j=0;j<param;j++){
			param_x_vec(j) = xvec_test_mat(i,j);
		}
		yvec_test_mat(i,0) = function(param_x_vec);
	}

	//Execute posterior function generation
	output_mat = gaussian_process_solver_regression(kernal_square_exp_noise_function, xvec_train_mat, yvec_train_mat, xvec_test_mat, samples, hyperp_vec, epsilon, regression_polynomial_function, power);

	//Write the output file
	write_output(output_mat, test, param, observables, samples, "posterior.dat");

	//Write to file as:
	//		param1, param2, ... , training value
	ofile.open("trainset.dat");
	for(i=0;i<train;i++){
		for(j=0;j<param;j++){
			ofile << std::setw(15) << xvec_train_mat(i,j);
		}
		ofile << std::setw(15) << yvec_train_mat(i,0) << std::endl;
	}
	ofile.close();


/*
	int resolution = 50,
		selection = 2,
		max_index1, max_index2;

	double theta_i = 0.8,
			theta_f = 0.9,
			dtheta = (theta_f - theta_i)/(double) (resolution - 1),
			precision = 0.0001;

	arma::mat hyperp_range_vec = arma::zeros<arma::mat> (2, hyperp);

	arma::vec  log_likelihood_vec = arma::zeros<arma::mat> (resolution),
				log_likelihood_der_vec = log_likelihood_vec;

	hyperp_range_vec(0, selection) = theta_i;
	hyperp_range_vec(1, selection) = theta_f;

	log_likelihood_vec = log_likelihood_function(yvec_train_mat, xvec_train_mat, resolution, kernal_square_exp_function, hyperp_vec, hyperp_range_vec, selection);
	max_index1 = find_max(log_likelihood_vec);
	log_likelihood_der_vec = log_likelihood_derivative(yvec_train_mat, xvec_train_mat, resolution, kernal_square_exp_function, hyperp_vec, hyperp_range_vec, selection, precision);
	max_index2 = find_zero(log_likelihood_der_vec);

	std::cout << std::setw(15) << theta_i + dtheta*(double) max_index1  << std::setw(15) << log_likelihood_vec(i) << std::endl;
	std::cout << std::setw(15) << theta_i + dtheta*(double) max_index2  << std::setw(15) << log_likelihood_der_vec(i) << std::endl;

	ofile.open("hyperp.dat");
	for(i=0;i<resolution;i++){
		for(j=0;j<1;j++){
			ofile << std::setw(15) << theta_i + dtheta*(double) i  << std::setw(15) << log_likelihood_vec(i) << std::setw(15) << log_likelihood_der_vec(i) << std::endl;
		}
	}
	ofile.close();
*/
	return 0;
}

