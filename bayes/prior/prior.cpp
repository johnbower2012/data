#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include<fstream>
#include<string>
#include "time.h"
#include "memory.cpp"
#include "armadillo"

std::ofstream ofile;

int main(int argc, char* argv[]){
	int i, j, k,
		length, input_values, samples;
	unsigned seed;
	double x_ini, x_fin, dx,
			xi, xj, sigma;

	std::string outfilename;

	clock_t start, finish;

	if(argc<8){
		std::cout << "Improper entry. Please also enter 'length test x_ini x_fin sigma samples outfilename' on same line." << std::endl;
		exit(1);
	}
	else{
		length = atoi(argv[1]);
		if(length<2){
			std::cout << "Improper entry. 'Length' must be greater than 1." << std::endl;
			exit(1);
		}
		else{
			length = atoi(argv[1]);
			test = atoi(argv[2]);
			x_ini = atof(argv[3]);
			x_fin = atof(argv[4]);
			sigma = atof(argv[5]);
			samples = atof(argv[6]);
			outfilename = argv[7];

			dx = (x_fin - x_ini)/(double) (length-1);
			input_values = 1;
			seed = std::chrono::system_clock::now().time_since_epoch().count();
		}
	}

	std::default_random_engine generator(seed);
	std::normal_distribution<double> dist(0,1.0);

	arma::mat 	xvec_train_mat = arma::zeros<arma::mat>(length,input_values),
				xvec_test_mat = arma::zeros<arma::mat>(length,input_values),
				cov_mat = arma::zeros<arma::mat>(length,length),
				L = cov_mat,
				
				prior_func = arma::zeros<arma::mat>(length,samples),
				random_sample = prior_func,

				I = arma::eye<arma::mat>(length,length);


	for(i=0;i<length;i++){
		for(j=0;j<input_values;j++){
			xvec_train_mat(i,j) = x_ini + ((double) i)*dx;
		}
	}

	dx = (x_fin - x_ini)/(double) (test+1);
	for(i=0;i<test;i++){
		for(j=0;j<input_values;j++){
			xvec_test_mat(i,j) = x_ini + ((double) (i+1))*dx;
		}
	}


	//calculate covariance matrix
	for(i=0;i<length;i++){
		xi = xvec_mat(i,0);
		for(j=0;j<length;j++){
			xj = xvec_mat(j,0);			
			cov_mat(i,j) = exp(-(xi-xj)*(xi-xj)/(2.0*sigma*sigma));
		}
	}

	//Cholesky Decomposition
	L = arma::chol(cov_mat + I*1e-15, "lower");

	//prior distribution
	for(i=0;i<length;i++){
		for(j=0;j<samples;j++){
			random_sample(i,j) = dist(generator);
		}
	}

	//sample prior functions
	prior_func = L*random_sample;

	//write to file
	ofile.open(outfilename);
	for(i=0;i<length;i++){
		ofile << std::setw(15) << xvec_mat(i,0);
		for(j=0;j<samples;j++){
			ofile << std::setw(15) << prior_func(i,j);
		}
		ofile << std::endl;
	}
	ofile.close();

	return 0;
}
