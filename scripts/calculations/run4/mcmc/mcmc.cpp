#include<iostream>
#include<iomanip>
#include<cmath>
#include<random>
#include<chrono>
#include<armadillo>

class MCMC{
public:
  int parameter_count;
  arma::vec position;
  arma::vec test_position;
  arma::vec widths;

  unsigned seed;
  std::default_random_engine generator;

  MCMC();
  void set_parameter_count(int new_parameter_count);
  void set_position(arma::vec new_position);
  void set_widths(arma::vec new_widths);

  void set_seed(int seed);
  void set_seed_clock();
  double normal();
  double uniform();

  void step();
  double get_likelihood();
};

int main(int argc, char* argv[]){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> dist(0,1.0);

  int parameters = 2;
  arma::vec position = arma::zeros<arma::vec>(parameters);
  arma::vec widths = arma::zeros<arma::vec>(parameters);
  for(int i=0;i<parameters;i++){
    position(i) = i+2;
    widths(i) = (double) (i+2)/(2.0);
  }
  MCMC random;
  random.set_parameter_count(parameters);
  random.set_position(position);
  random.set_widths(widths);

  random.position.print("p");
  for(int i=0;i<1000;i++){
    random.step();
  }
  random.position.print("p");  
  return 0;
}

/*SETUP--
***************/
MCMC::MCMC(){
  set_seed(1);
}
void MCMC::set_parameter_count(int new_parameter_count){
  parameter_count = new_parameter_count;
  position = arma::zeros<arma::vec>(parameter_count);
  widths = arma::zeros<arma::vec>(parameter_count);
  test_position = arma::zeros<arma::vec>(parameter_count);
}
void MCMC::set_position(arma::vec new_position){
  position = new_position;
}
void MCMC::set_widths(arma::vec new_widths){
  widths = new_widths;
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
void MCMC::step(){
  double random, LH;
  for(int i=0;i<parameter_count;i++){
    test_position(i) = position(i) + normal()*widths(i);
  }
  LH = get_likelihood();

  if(LH<1.0){
    random = uniform();
    if(random<LH){
      position = test_position;
    }
  } else{
    position = test_position;
  }
}      
