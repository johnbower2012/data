#include "pca.hpp"

//calculate the means of an assembly of column-wise vectors
arma::vec calculate_mean_function(arma::mat val_matrix){
	int i, j, 
		rows = val_matrix.n_rows,
		cols = val_matrix.n_cols;
	arma::vec mean_vec = arma::zeros<arma::vec> (cols);

	for(i=0;i<cols;i++){
		for(j=0;j<rows;j++){
			mean_vec(i) += val_matrix(j,i);
		}
		mean_vec(i) /= (double) rows;
	}
	return mean_vec;
}
//creates covariance matrix for a matrix of n-repitions of p-dim row vectors
arma::mat calculate_covariance_function(arma::mat val_matrix, arma::vec mean_vec){
	int i, j, k,
		reps = val_matrix.n_rows,
		dims = val_matrix.n_cols;
	double xi_mean, xj_mean, x_ki, x_kj;
	arma::mat cov_matrix = arma::zeros<arma::mat> (dims,dims);

	for(i=0;i<dims;i++){
		xi_mean = mean_vec(i);
		for(j=0;j<dims;j++){
			xj_mean = mean_vec(j);
			for(k=0;k<reps;k++){
				x_ki = val_matrix(k,i);
				x_kj = val_matrix(k,j);
				cov_matrix(i,j) += (x_ki - xi_mean)*(x_kj - xj_mean);
			}
			cov_matrix(i,j) /= (double) (reps);
		}
	}
	return cov_matrix;
}
//calculate zmatrix, with mean zero and unit std, i.e. [z = (y-ybar)/ystd]
arma::mat calculate_zmatrix_function(arma::mat val_matrix, arma::vec mean_vec, arma::mat cov_matrix){
	int i, j,
		reps = val_matrix.n_rows,
		dims = val_matrix.n_cols;
	double xi_mean, xi_std;
	arma::mat zval_matrix = arma::zeros<arma::mat> (reps,dims);

	for(i=0;i<dims;i++){
		xi_mean = mean_vec(i);
		xi_std = sqrt(cov_matrix(i,i));
		for(j=0;j<reps;j++){
			zval_matrix(j,i) = (val_matrix(j,i) - xi_mean)/xi_std;
		}
	}
	return zval_matrix;
}
//Calculate valmatrix from zmatrix
arma::mat calculate_valmatrix_function(arma::mat zval_matrix, arma::vec mean_vec, arma::mat cov_matrix){
	int i, j,
		reps = zval_matrix.n_rows,
		dims = zval_matrix.n_cols;
	double xi_mean, xi_std;
	arma::mat val_matrix = arma::zeros<arma::mat> (reps,dims);

	for(i=0;i<dims;i++){
		xi_mean = mean_vec(i);
		xi_std = sqrt(cov_matrix(i,i));
		for(j=0;j<reps;j++){
			val_matrix(j,i) = (zval_matrix(j,i)*xi_std + xi_mean);
		}
	}
	return val_matrix;
}
//creates zcovariance matrix for a matrix of n-repitions of p-dim row vectors, with zero column-wise mean
arma::mat calculate_zcovariance_function(arma::mat zval_matrix){
	int i, j, k,
		reps = zval_matrix.n_rows,
		dims = zval_matrix.n_cols;
	double x_ki, x_kj;
	arma::mat zcov_matrix = arma::zeros<arma::mat> (dims,dims);

	for(i=0;i<dims;i++){
		for(j=0;j<dims;j++){
			for(k=0;k<reps;k++){
				x_ki = zval_matrix(k,i);
				x_kj = zval_matrix(k,j);
				zcov_matrix(i,j) += x_ki*x_kj;
			}
			zcov_matrix(i,j) /= (double) (reps);
		}
	}
	return zcov_matrix;
}
//sort eigenvectors according to eigenvalues, biggest to smallest --> L1 > L2 > ... > Ln
arma::mat sort_eigenvectors_function(arma::vec eigval_vec, arma::mat eigvec_matrix){
	int i, j, sort,
		dims = eigval_vec.n_elem;
	arma::mat 	eigsort_matrix = arma::zeros<arma::mat> (dims,dims),
				eigval_matrix = arma::zeros<arma::mat> (dims,2);

	for(i=0;i<dims;i++){
		eigval_matrix(i,0) = eigval_vec(i);
		eigval_matrix(i,1) = i;
	}
	eigval_matrix = arma::sort(eigval_matrix, "descending", 0);

	for(i=0;i<dims;i++){
		sort = eigval_matrix(i,1);
		for(j=0;j<dims;j++){
			eigsort_matrix(j,i) = eigvec_matrix(j,sort);
		}
	}
	return eigsort_matrix;
}
//full transformation from val_matrix to tval_matrix
arma::mat calculate_tmatrix_function(arma::mat val_matrix){
	int reps = val_matrix.n_rows,
		dims = val_matrix.n_cols;

	arma::mat	zval_matrix = arma::zeros<arma::mat>(reps,dims),
				tval_matrix = zval_matrix,

				cov_matrix = arma::zeros<arma::mat>(dims,dims),
				zcov_matrix = cov_matrix,
				eigvec_matrix = cov_matrix,
				eigsort_matrix = cov_matrix,

				eigval_matrix = arma::zeros<arma::mat>(dims,2);

	arma::vec	mean_vec = arma::zeros<arma::vec>(dims),
				eigval_vec = mean_vec;

	mean_vec = calculate_mean_function(val_matrix);
	cov_matrix = calculate_covariance_function(val_matrix, mean_vec);
	zval_matrix = calculate_zmatrix_function(val_matrix, mean_vec, cov_matrix);
	zcov_matrix = calculate_zcovariance_function(zval_matrix);

	//find eigenvalues and eigenvectors for zcov matrix
	arma::eig_sym(eigval_vec, eigvec_matrix, zcov_matrix);

	//sort eigenvalues in descending order
	eigvec_matrix = sort_eigenvectors_function(eigval_vec, eigvec_matrix);

	//calculate rotated system
	tval_matrix = zval_matrix*eigvec_matrix;

	return tval_matrix;
}
//full transformation from val_matrix to tval_matrix, return mean & cov matrices
arma::mat calculate_tmatrix_function(arma::mat val_matrix, arma::mat& eigvec_matrix, arma::mat& mean_vec, arma::mat& cov_matrix){
	int reps = val_matrix.n_rows,
		dims = val_matrix.n_cols;

	arma::mat	zval_matrix = arma::zeros<arma::mat>(reps,dims),
				tval_matrix = zval_matrix,

				zcov_matrix = cov_matrix,
				eigsort_matrix = cov_matrix,

				eigval_matrix = arma::zeros<arma::mat>(dims,2);

	arma::vec	eigval_vec = arma::zeros<arma::vec>(dims);

	cov_matrix = arma::zeros<arma::mat>(dims,dims);
	eigvec_matrix = arma::zeros<arma::mat>(dims,dims);
	mean_vec = arma::zeros<arma::vec>(dims);

	mean_vec = calculate_mean_function(val_matrix);
	cov_matrix = calculate_covariance_function(val_matrix, mean_vec);
	zval_matrix = calculate_zmatrix_function(val_matrix, mean_vec, cov_matrix);
	zcov_matrix = calculate_zcovariance_function(zval_matrix);

	//find eigenvalues and eigenvectors for zcov matrix
	arma::eig_sym(eigval_vec, eigvec_matrix, zcov_matrix);

	//sort eigenvalues in descending order
	eigvec_matrix = sort_eigenvectors_function(eigval_vec, eigvec_matrix);

	//calculate rotated system
	tval_matrix = zval_matrix*eigvec_matrix;

	return tval_matrix;
}
//Transformation from zval_matrix to tval_matrix
arma::mat calculate_ztotmatrix_function(arma::mat zval_matrix){
	int reps = zval_matrix.n_rows,
		dims = zval_matrix.n_cols;

	arma::mat	tval_matrix = arma::zeros<arma::mat>(reps,dims),

				cov_matrix = arma::zeros<arma::mat>(dims,dims),
				zcov_matrix = cov_matrix,
				eigvec_matrix = cov_matrix,
				eigsort_matrix = cov_matrix,

				eigval_matrix = arma::zeros<arma::mat>(dims,2);

	arma::vec	mean_vec = arma::zeros<arma::vec>(dims),
				zmean_vec = mean_vec,
				eigval_vec = mean_vec;


	zcov_matrix = calculate_zcovariance_function(zval_matrix);

	//find eigenvalues and eigenvectors for zcov matrix
	arma::eig_sym(eigval_vec, eigvec_matrix, zcov_matrix);

	//sort eigenvalues in descending order
	eigvec_matrix = sort_eigenvectors_function(eigval_vec, eigvec_matrix);

	//calculate rotated system
	tval_matrix = zval_matrix*eigvec_matrix;

	return tval_matrix;
}
//full transformation from tval_matrix to val_matrix
arma::mat calculate_ttovalmatrix_function(arma::mat tval_matrix, arma::mat eigvec_matrix, arma::vec mean_vec, arma::mat cov_matrix){
	int i, j,
		reps = tval_matrix.n_rows,
		dims = tval_matrix.n_cols;

	double xi_mean, xi_std;

	arma::mat	val_matrix = arma::zeros<arma::mat> (reps,dims);

	//calculate rotated system
	val_matrix = tval_matrix*eigvec_matrix.t();

	//calculate system including old mean and non-unit std
	for(i=0;i<dims;i++){
		xi_mean = mean_vec(i);
		xi_std = sqrt(cov_matrix(i,i));
		for(j=0;j<reps;j++){
			val_matrix(j,i) = val_matrix(j,i)*xi_std + xi_mean;
		}
	}

	return val_matrix;
}


/*****Writing to File
******************************************/
void writeFile(std::string filename, arma::mat val_matrix){
	std::fstream ofile;
	int reps = val_matrix.n_rows,
		dims = val_matrix.n_cols;
	ofile.open(filename);
	for(int i=0;i<reps;i++){
		for(int j=0;j<dims;j++){
			ofile << std::setw(15) << val_matrix(i,j);
		}
		ofile << std::endl;
	}
	ofile.close();
}




