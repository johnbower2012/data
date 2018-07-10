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
#include "pca.cpp"
#include "gaussianprocess.cpp"

std::ofstream ofile;

int main(int argc, char* argv[]){
	int 	reps, dims,
			i, j;
	unsigned seed;
	double random;

	if(argc<3){
		std::cout << "Improper input. Enter also 'dimensionality repetitions' on same line." << std::endl;
		exit(1);
	}
	else{
		dims = atoi(argv[1]);
		reps = atoi(argv[2]);
		seed = std::chrono::system_clock::now().time_since_epoch().count();
	}

//random numbers
	std::default_random_engine generator (seed);
	std::normal_distribution<double> dist(0.0,0.5);

//generate matrices and arrays
	arma::mat	val_matrix = arma::zeros<arma::mat>(reps,dims),
				zval_matrix = val_matrix,
				tval_matrix = val_matrix,

				cov_matrix = arma::zeros<arma::mat>(dims,dims),
				zcov_matrix = cov_matrix,
				eigvec_matrix = cov_matrix;

	arma::vec	mean_vec = arma::zeros<arma::mat>(dims),
				eigval_vec = mean_vec;

//initialize matrices and arrays
	for(i=0;i<reps;i++){
		for(j=0;j<dims;j++){
			random = dist(generator);
			val_matrix(i,j) = (double) ((i+1)*(j+1))*0.1 + random;
			if(j==dims-1){
				val_matrix(i,j) = random;
			}
		}
	}
	

//calculate the PCA Transformed matrix
	mean_vec = calculate_mean_function(val_matrix);
	cov_matrix = calculate_covariance_function(val_matrix, mean_vec);
	zval_matrix = calculate_zmatrix_function(val_matrix, mean_vec, cov_matrix);

	zcov_matrix = calculate_zcovariance_function(zval_matrix);
	arma::eig_sym(eigval_vec, eigvec_matrix, zcov_matrix);
	eigvec_matrix = sort_eigenvectors_function(eigval_vec, eigvec_matrix);
	eigval_vec.print();
	eigvec_matrix.print();

	tval_matrix = zval_matrix*eigvec_matrix;

	double sum;
	int k;
	for(i=0;i<dims;i++){
		for(j=0;j<dims;j++){
			sum = 0.0;
			for(k=0;k<reps;k++){
				sum += tval_matrix(k,i)*tval_matrix(k,j);
			}
			sum /= (double) reps;
			std::cout << std::setw(10) << i << std::setw(10) << j << std::setw(15) << sum << std::endl;
		}
	}
/******
 *
 *WRITE TO FILE
 *
 ******/

//original values
	ofile.open("val.dat");
	for(i=0;i<reps;i++){
		for(j=0;j<dims;j++){
			ofile << std::setw(15) << val_matrix(i,j);
		}
		ofile << std::endl;
	}
	ofile.close();
//PCA transformed z-system, along principle components
	ofile.open("tval.dat");
	for(i=0;i<reps;i++){
		for(j=0;j<dims;j++){
			ofile << std::setw(15) << tval_matrix(i,j);
		}
		ofile << std::endl;
	}
	ofile.close();


	return 0;

}
