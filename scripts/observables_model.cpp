#include<iostream>
#include<fstream>
#include<cmath>
#include "armadillo"
#include "pricomana.cpp"
#include "observables_class.cpp"
void load_file(std::string fileName, arma::mat &file){
  int load = file.n_rows;
  int param = file.n_cols;
  std::ifstream ifile;
  ifile.open(fileName);
  for(int i=0;i<load;i++){
    for(int j=0;j<param;j++){
      ifile >> file(i,j);
    }
  }
  ifile.close();
}
int main(int argc, char* argv[]){
  std::string *infilename,*outfilename;
  int observables,obs_file,files,lines,runs;
  int i,j,k;
  infilename = new std::string[observables];
  outfilename = new std::string[observables];
  std::string dest_folder;

  files=4;
  obs_file=3;
  observables=files*obs_file;

  if(argc<4+files){
    std::cout << "Improper input. Enter also 'lines runs observables ifn*[7]' on same line." << std::endl;
    exit(1);
  }
  else{
    infilename = new std::string[files];
    dest_folder = argv[1];
    lines=atoi(argv[2]);
    runs=atoi(argv[3]);
    files=atoi(argv[4]);
    for(i=0;i<files;i++){
      infilename[i]=argv[5+i];
    }
  }
  printf("Arguments read in as:\n");
  printf("observables.x %s %d %d %d",dest_folder.c_str(),lines,runs,files);
  for(i=0;i<files;i++){
    printf(" %s",infilename[i].c_str());
  }
  printf("\n");

  /*********
	     LOAD FILE
	     CALCULATE OBS_MAT
  *********/
  arma::mat *val_matrix, obs_matrix;
  arma::vec delY_vec, obs_error;
  val_matrix = new arma::mat[files];
  obs_matrix = arma::zeros<arma::mat>(runs,observables);
  delY_vec = arma::zeros<arma::vec>(lines);
  obs_error = arma::zeros<arma::vec>(observables);
  //all observables that aren't pipi, reduce unc by factor of ten.

  /*
  //pipi
  int part=0;
  obs_error(part*obs_file) = 0.01;
  obs_error(1+part*obs_file) = 0.01;
  obs_error(2+part*obs_file) = 0.005;
  //ppbar
  part=1;
  obs_error(part*obs_file) = 0.005;
  obs_error(1+part*obs_file) = 0.005;
  obs_error(2+part*obs_file) = 0.001;
  //pK
  part=2;
  obs_error(part*obs_file) = 0.001;
  obs_error(1+part*obs_file) = 0.001;
  obs_error(2+part*obs_file) = 0.0005;
  //KK
  part=3;
  obs_error(part*obs_file) = 0.005;
  obs_error(1+part*obs_file) = 0.005;
  obs_error(2+part*obs_file) = 0.001;
  */

  load_file(files, lines, runs, infilename, delY_vec, val_matrix);
  obs_matrix_moments_fabs(files,obs_file,val_matrix,delY_vec,obs_matrix);

  /*********
	     CONDUCT PCA
	     PRINT
  *********/

  arma::mat tval_matrix,zcov_matrix,eigvec_matrix,print_matrix;
  arma::vec eigval_vec,mean_vec;
  print_matrix = arma::zeros<arma::mat>(observables,observables+1);
  tval_matrix = arma::zeros<arma::mat>(observables,observables);
  zcov_matrix = arma::zeros<arma::mat>(observables,observables);
  eigvec_matrix = arma::zeros<arma::mat>(observables,observables);
  eigval_vec = arma::zeros<arma::vec>(observables);
  mean_vec = arma::zeros<arma::vec>(observables);

  arma::mat err = arma::zeros<arma::mat>(1,observables);
  arma::mat y_tilde = arma::zeros<arma::mat>(lines,observables);
  load_file("moments_exp_data.dat",err);
  for(int i=0;i<observables;i++){
    obs_error(i) = 0.1*err(0,i);
  }
  tilde_function(obs_matrix,obs_error,mean_vec,y_tilde);
  tval_matrix = calculate_tmatrix_function(obs_matrix, obs_error, eigval_vec, eigvec_matrix, mean_vec, zcov_matrix);

  print_fractional_sum(eigval_vec);
  eigval_vec.print("+++++ eigvalues +++++");
  eigvec_matrix.print("+++++ eigvectors +++++");

  for(i=0;i<observables;i++){
    print_matrix(i,0) = eigval_vec(i);
    for(j=0;j<observables;j++){
      print_matrix(i,j+1) = eigvec_matrix(j,i);
    }
  }
  std::string printname;
  std::string title;
  printname = "moments_model_data.dat";
  title = "#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,obs_matrix);

  printname = "moments_model_data_tilde.dat";
  title = "#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,y_tilde);

  printname = "moments_errors.dat";
  title="#pipi/ppbar/pK/KK";
  print_file(printname,obs_error);

  printname = "moments_means.dat";
  title="#pipi/ppbar/pK/KK";
  print_file(printname,mean_vec);

  printname = "moments_model_eigvec.dat";
  title = "#colvec\n#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,eigvec_matrix);

  printname = "moments_model_pca.dat";
  title = "#eigval col -- eigvec rows\n#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,print_matrix);

  print_matrix = tval_matrix;
  printname = "moments_model_z.dat";
  title = "#y~*eigvec";
  print_file(printname,print_matrix);

  delete[] val_matrix;
  delete[] infilename;
  delete[] outfilename;

  return 0;
}

