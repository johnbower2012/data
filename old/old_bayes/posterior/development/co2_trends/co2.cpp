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
std::ifstream ifile;

double function(arma::vec param){
	double value = param(0) + param(0)*sin(1.5*param(0));
	return value;
}
bool identical(arma::vec test){
	int i, j, length = test.n_elem;
	bool value = false;
	for(i=0;i<length;i++){
		for(j=i+1;j<length;j++){
			if(test(i)==test(j)){
				value = true;
				break;
			}
		}
	}
	return value;
}
arma::vec hyperp_optimize_function(
	arma::mat Y_mat, arma::mat X_mat, 
	int resolution, double precision, 
	arma::mat kernal_func(arma::mat, arma::mat, arma::vec), 
	arma::vec hyperp_guess_vec, arma::mat& hyperp_range_mat);

int main(int argc, char* argv[]){
	//Variables for use throughout the program
	int i,
		train, test, 
		param, observables, hyperp,
		samples;
	unsigned seed;
	double 	x, epsilon;

	//For writing to file
	std::string infilename, outfilename;

	//Test to ensure the proper input. Terminate program if command line inputs are not given
	if(argc<3){
		std::cout << "Improper entry. Please also enter 'infilename outfilename' on same line." << std::endl;
		exit(1);
	}
	else{
		infilename = argv[1];
		outfilename = argv[2];
	}

	/*SYSTEM and OUTPUT--
	***********************************************************/

		train = 550;
		test = 20*12;

		samples = 3;
		param = 1;
		observables = 1;

		epsilon = 1e-8;

		seed = std::chrono::system_clock::now().time_since_epoch().count();

	/*HYPERPARAMETERS--
	***********************************************************/

		hyperp = 11;

	//random number generator
	std::default_random_engine generator(seed);
	std::normal_distribution<double> dist(0,1.0);

	//armadillo matrices and vectors for computational use
	arma::mat 	xvec_train_mat = arma::zeros<arma::mat>(train,param),
				yvec_train_mat = arma::zeros<arma::mat>(train,observables),
				xvec_test_mat = arma::zeros<arma::mat>(test,param),

				output_mat = arma::zeros<arma::mat>(test,param + 3 + samples);

	arma::vec	hyperp_vec = arma::zeros<arma::vec> (hyperp),
				param_x_vec = arma::zeros<arma::vec> (param),
				years_vec = arma::zeros<arma::vec> (train),
				months_vec = years_vec,
				co2_vec = years_vec;

	//hyperparamter vector for use in 'function'

	hyperp_vec(0) = 66.0;
	hyperp_vec(1) = 67.0;
	hyperp_vec(2) = 2.4;
	hyperp_vec(3) = 90.0;
	hyperp_vec(4) = 1.3;
	hyperp_vec(5) = 0.66;
	hyperp_vec(6) = 1.2;
	hyperp_vec(7) = 0.78;
	hyperp_vec(8) = 0.18;
	hyperp_vec(9) = 1.6/12.0;
	hyperp_vec(10) = 0.19;
/*
	hyperp_vec(0) = 250.24;
	hyperp_vec(1) = 64.72;
	hyperp_vec(2) = 2.08;
	hyperp_vec(3) = 83.2;
	hyperp_vec(4) = 1.4;
	hyperp_vec(5) = 0.72;
	hyperp_vec(6) = 1.216;
	hyperp_vec(7) = 0.488;
	hyperp_vec(8) = 0.184;
	hyperp_vec(9) = 0.128;
	hyperp_vec(10) = 0.1936;
*/
	//Generate Functions	
	ifile.open(infilename);
	for(i=0;i<train;i++){
		ifile >> x; ifile >> months_vec(i); ifile >> years_vec(i); ifile >> x; ifile >> co2_vec(i);
		ifile >> x; ifile >> x;
		xvec_train_mat(i,0) = years_vec(i);
		yvec_train_mat(i,0) = co2_vec(i);
	}
	ifile.close();
	for(i=0;i<test;i++){
		xvec_test_mat(i,0) = xvec_train_mat(train-1,0) + (double) (i)/12.0;
	}

	//Execute posterior function generation
	output_mat = gaussian_process_solver_basic(kernal_co2_function, xvec_train_mat, yvec_train_mat, xvec_test_mat, samples, hyperp_vec, epsilon);

	//Write the output file
	write_output(output_mat, test, param, observables, samples, outfilename);

	//Write the trainingset to file
	write_trainset(xvec_train_mat, yvec_train_mat, "trainset.dat");

/*
	int resolution = 11;
	double precision = 0.0001,
			theta_i = 1e-12,
			theta_f = 10.0;

	arma::mat hyperp_range_mat = arma::zeros<arma::mat>(2,hyperp);
	arma::vec hyperp_guess_vec = arma::zeros<arma::mat>(hyperp);

	for(i=0;i<hyperp;i++){
		hyperp_range_mat(0,i) = theta_i;
		hyperp_range_mat(1,i) = theta_f;
	}

	hyperp_guess_vec(0) = 250.24;
	hyperp_guess_vec(1) = 64.72;
	hyperp_guess_vec(2) = 2.08;
	hyperp_guess_vec(3) = 83.2;
	hyperp_guess_vec(4) = 1.4;
	hyperp_guess_vec(5) = 0.72;
	hyperp_guess_vec(6) = 1.216;
	hyperp_guess_vec(7) = 0.488;
	hyperp_guess_vec(8) = 0.184;
	hyperp_guess_vec(9) = 0.128;
	hyperp_guess_vec(10) = 0.1936;

hyperp_vec = hyperp_optimize_function(yvec_train_mat, xvec_train_mat, resolution, precision, kernal_co2_function, hyperp_guess_vec, hyperp_range_mat);

hyperp_vec.print();
*/
	return 0;
}

arma::vec hyperp_optimize_function(
	arma::mat Y_mat, arma::mat X_mat, 
	int resolution, double precision, 
	arma::mat kernal_func(arma::mat, arma::mat, arma::vec), 
	arma::vec hyperp_guess_vec, arma::mat& hyperp_range_mat){
	
	int i,
		zeros, selection,
		hyperp = hyperp_guess_vec.n_elem,
		count = 0;

	double theta_i, theta_f, dtheta, scale;

	arma::vec ll_der_vec = arma::zeros<arma::vec>(resolution),
				zeros_vec,
				hyperp_values_vec,
				ll_vec,
				hyperp_vec = hyperp_guess_vec,
				temp_vec;

	for(i=0;i<hyperp;i++){
			hyperp_vec(i) = hyperp_guess_vec(i);
	}

	for(i=0;i<hyperp;i++){
		selection = i;
		zeros=0;
		std::cout << "hyperparameter" << selection << ":" << std::endl;
		count=0;
		
	hyperp_vec.print();
		while(zeros==0&count<5){
			ll_der_vec =  log_likelihood_derivative(Y_mat, X_mat, resolution, kernal_func, hyperp_vec, hyperp_range_mat, selection, precision);
			zeros_vec = find_zeros(ll_der_vec);
			zeros = zeros_vec.n_elem-1;
			if(zeros==1){
				theta_i = hyperp_range_mat(0,selection);
				theta_f = hyperp_range_mat(1,selection);
				dtheta = (theta_f - theta_i)/(double) (resolution-1);
				scale = (theta_f - theta_i)/theta_f;
				count = 0;
				while(scale>0.01){
					theta_i = hyperp_range_mat(0,selection);
					theta_f = hyperp_range_mat(1,selection);
					dtheta = (theta_f - theta_i)/(double) (resolution-1);

					std::cout << "hyperp" << selection << ":" << std::endl;
					std::cout << std::setw(15) << "zeros: " << zeros_vec(0) << std::setw(15) << "hyperp:" << hyperp_values_vec(0) << std::endl;					
					std::cout << std::setw(15) << "searching in " << theta_i << " to " << theta_f << std::endl;
					
					ll_der_vec =  log_likelihood_derivative(Y_mat, X_mat, resolution, kernal_func, hyperp_vec, hyperp_range_mat, selection, precision);
					zeros_vec = find_zeros(ll_der_vec);
					zeros = zeros_vec.n_elem-1;

					theta_i += dtheta*(double) (zeros_vec(0) - 1);
					theta_f = theta_i + 2.0*dtheta;

					hyperp_range_mat(0,selection) = theta_i;
					hyperp_range_mat(1,selection) = theta_f;
					scale = (theta_f - theta_i)/theta_f;

					std::cout << std::setw(15) << "new zero at " << theta_i + dtheta;
				}
				hyperp_vec(selection) = theta_i + dtheta;
			}
			else if(zeros>1){
				std::cout << "Multiple zeros found for hyperparamter" << selection << ", skipping to next." << std::endl;
			}
			else{
				hyperp_range_mat(1,selection) *= 2.0;
				count++;
				if(count==5){
					std::cout << "Search range reached, skipping to next hyperparameter." << std::endl;
				}
			}
		}
	}
	return hyperp_vec;
}







