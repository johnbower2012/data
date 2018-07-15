#include<iostream>
#include<iomanip>
#include<cmath>
#include<random>
#include<chrono>
#include<armadillo>
#include<string>
#include "emulator_class.cpp"

class MCMC{
public:
  double target_value;
  int parameter_count;
  double sigma;
  arma::vec position;
  arma::vec test_position;
  arma::vec widths;

  unsigned seed;
  std::default_random_engine generator;

  MCMC();
  void set_target_value(double target_value_new);
  void set_sigma(double new_sigma);
  void set_parameter_count(int new_parameter_count);
  void set_position(arma::vec new_position);
  void set_widths(arma::vec new_widths);

  void set_seed(int seed);
  void set_seed_clock();
  double normal();
  double uniform();

  void step();
  double get_likelihood();
  double get_likelihood(double z);
  bool decide(double z);
};

int main(int argc, char* argv[]){
  int mcmc_runs = 100;
  std::string datafile, expfile, hpfile, betafile, outfile;

  if(argc<7){
    printf("Improper usage. Please enter also 'model_file exp_file hpfile betafile mcmc_runs outfile' on sabme line.\n");
    exit(1);
  }else{
    datafile = argv[1];
    expfile = argv[2];
    hpfile = argv[3];
    betafile = argv[4];
    mcmc_runs = atoi(argv[5]);
    outfile = argv[6];
  }
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::normal_distribution<double> dist(0.0,1.0);
  std::default_random_engine generator(seed);

  double target_value;
  double epsilon=1e-8;
  int parameters = 4;
  int observables = 12;
  int model_runs = 1000;
  int exp_runs = 1;
  int hyperparameteres = 3;
  arma::mat X = arma::zeros<arma::mat>(model_runs,parameters);
  arma::mat Y_model = arma::zeros<arma::mat>(model_runs,observables);
  arma::mat y_model = arma::zeros<arma::mat>(model_runs,1);
  arma::mat Y_exp = arma::zeros<arma::mat>(exp_runs,observables);
  arma::vec hyp = arma::zeros<arma::vec>(hyperparameteres);
  arma::mat beta = arma::zeros<arma::mat>(observables,parameters+1);
  arma::mat X_s = arma::zeros<arma::mat>(1,parameters);
  load_data_file(datafile,X,Y_model);
  load_file(expfile,Y_exp);
  load_file(hpfile,hyp);
  load_beta_file(betafile,beta);
  y_model = Y_model.col(0);
  target_value = Y_exp(0,0);

  arma::vec position = arma::zeros<arma::vec>(parameters);
  arma::vec widths = arma::zeros<arma::vec>(parameters);
  for(int i=0;i<parameters;i++){
    position(i) = dist(generator);
    widths(i) = (1.5-0.1);
  }

  emulator gauss(X,hyp,beta,epsilon);
  MCMC random;
  random.set_target_value(target_value);
  random.set_sigma(0.001);
  random.set_parameter_count(parameters);
  random.set_position(position);
  random.set_widths(widths);
  random.set_seed(1);

  arma::mat gauss_mat = arma::zeros<arma::mat>(1,parameters+3);
  double z=0.0;
  bool stepped;
  double fraction=0.0;
  int trace_count = (int) (0.2*mcmc_runs);
  arma::mat trace = arma::zeros<arma::mat>(trace_count,parameters);
  int j=0;
  for(int i=0;i<mcmc_runs;i++){
    random.step();
    X_s = random.test_position.t();
    gauss_mat = gauss.emulate(X_s,y_model);
    z = gauss_mat(0,parameters);
    stepped = random.decide(z);
    printf("TV: %f Z: %f\n", random.target_value, z);
    if(stepped==true){
      fraction+=1.0;
    }
    if(i%5==0){
      trace.row(j) = random.position.t();
      j++;
    }
  }
  printf("fraction: %f\n",fraction/mcmc_runs);

  trace.print("trace");

  return 0;
}

/*SETUP--
***************/
MCMC::MCMC(){
  this->sigma = 1.0;
  set_seed(1);
}
void MCMC::set_target_value(double new_target_value){
  this->target_value = new_target_value;
}
void MCMC::set_sigma(double new_sigma){
  this->sigma = new_sigma;
}
void MCMC::set_parameter_count(int new_parameter_count){
  this->parameter_count = new_parameter_count;
  this->position = arma::zeros<arma::vec>(this->parameter_count);
  this->widths = arma::zeros<arma::vec>(this->parameter_count);
  this->test_position = arma::zeros<arma::vec>(this->parameter_count);
}
void MCMC::set_position(arma::vec new_position){
  this->position = new_position;
}
void MCMC::set_widths(arma::vec new_widths){
  this->widths = new_widths;
}

/*RANDOM--
***************/
void MCMC::set_seed(int seed){
  generator.seed(seed);
}
void MCMC::set_seed_clock(){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator.seed(seed);
}
double MCMC::normal(){
  std::normal_distribution<double> normal_dist(0.0,1.0);
  return normal_dist(generator);
}
double MCMC::uniform(){
  std::uniform_real_distribution<double> uniform_dist(0.0,1.0);
  return uniform_dist(generator);
}

/*STEP--
**************/
double MCMC::get_likelihood(){
  double z1 = 0.0, z2 = 0.0;
  for(int i=0;i<parameter_count;i++){
    z1 += position(i)*position(i);
    z2 += test_position(i)*test_position(i);
  }
  return exp(z1-z2);
}
double MCMC::get_likelihood(double z){
  double diff = z - target_value;
  return exp(-diff*diff/(2.0*sigma));
}

void MCMC::step(){
  double random, LH;
  for(int i=0;i<parameter_count;i++){
    test_position(i) = position(i) + normal()*widths(i);
  }
}      
bool MCMC::decide(double z){
  double random, LH;
  LH = get_likelihood(z);

  if(LH<1.0){
    random = uniform();
    if(random<LH){
      position = test_position;
      return true;
    }else{
      return false;
    }
  } else{
    position = test_position;
    return true;
  }
}      
/*FILEWRITING--
**************/

