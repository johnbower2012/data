#include<iostream>
#include<cmath>
#include<iomanip>
#include<random>
#include<armadillo>

/******************** 
  CLASS DECLARATION
 ********************/
class particle{
public:
  int dimensions;
  double mass;
  double charge;
  arma::vec position;
  arma::vec velocity;

  particle();
  particle(int Dimensions, double Mass, double Charge, arma::vec Position, arma::vec Velocity);
  
  void update(arma::vec position, arma::vec velocity);
  void print();
};

class particle_system{
public:
  int particle_count;
  int dimensions;

  double delta_t;
  double gamma;
  double sigma;

  std::vector<particle> system;
  arma::mat position;
  //  arma::vec (*force)(particle part1, particle part2);

  particle_system(int Dimensions, double Delta_t, double Gamma, double Sigma);
  
  void add(particle new_particle);
  arma::vec relative_coordinate();
  arma::vec langevin_force();
  arma::vec electrostatic_force();

  void verlet_step();

  void update();
  void run(int steps);
  void print();
};

/******************** 
  MAIN FUNCTION
 ********************/
int main(int argc, char* argv[]){
  double delta_t=0.1;
  double gamma=1.0;
  double sigma=1.0;

  int dimensions=2;
  int particles=3;
  particle_system system(dimensions, delta_t, gamma, sigma);
  arma::vec x = arma::zeros<arma::vec>(dimensions);
  arma::vec v = arma::zeros<arma::vec>(dimensions);
  arma::vec force_sum = arma::zeros<arma::vec>(dimensions);
  particle part1(dimensions,1,1,x,v);
  for(int i=0;i<particles;i++){
    part1.mass = 1.0;
    part1.charge = 1.0;
    for(int j=0;j<dimensions;j++){
      part1.position(j) = (i+j);
    }
    part1.velocity = v;
    system.add(part1);
  }

  system.position.print();
  arma::vec temp = system.relative_coordinate();
  temp.print();
  temp = system.electrostatic_force();
  temp.print();
  return 0;
}

/******************** 
  CLASS FUNCTIONS
 ********************/
/******************** 
  particle
 ********************/
particle::particle(){
  this->dimensions = 0;
  this->charge = 0;
  this->position = arma::zeros<arma::vec>(0);
  this->velocity = arma::zeros<arma::vec>(0);
}
particle::particle(int Dimensions, double Mass, double Charge, arma::vec Position, arma::vec Velocity){
  this->dimensions = Dimensions;
  this->mass = Mass;
  this->charge = Charge;
  this->position = Position;
  this->velocity = Velocity;
}
void particle::update(arma::vec Position, arma::vec Velocity){
  this->position = Position;
  this->velocity = Velocity;
}
void particle::print(){
  printf("dim: %d\ncharge: %f\n",this->dimensions,this->charge);
  this->position.print("position");
  this->velocity.print("velocity");
}

/******************** 
  particle_system
 ********************/

particle_system::particle_system(int Dimensions, double Delta_t, double Gamma, double Sigma){
  this->dimensions = Dimensions;
  this->delta_t = Delta_t;
  this->gamma = Gamma;
  this->sigma = Sigma;
  this->particle_count=0;
}
void particle_system::add(particle new_particle){
  this->system.push_back(new_particle);
  this->particle_count++;
  this->position.resize(dimensions,particle_count);
  this->position.col(particle_count-1) = new_particle.position;
}
arma::vec particle_system::relative_coordinate(){
  int index=0;
  int dim=this->dimensions+1;
  double r2=0.0;
  arma::vec relcoord = arma::zeros<arma::vec>(particle_count*(particle_count-1)*dim/2);
  for(int i=0;i<this->particle_count;i++){
    for(int j=0;j<i;j++){
      r2=0.0;
      for(int k=0;k<dim;k++){
	index = (i*(i-1)/2 + j)*(dim) + k;
	if(k==dim-1){
	  relcoord(index) = sqrt(r2);
	  continue;

	}
	else{
	  relcoord(index) = position(k,i) - position(k,j);
	  r2 += relcoord(index)*relcoord(index);
	}
      }
    }       
  }
  return relcoord;
}
arma::vec particle_system::langevin_force(){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  sigma = sqrt(delta_t);
  std::normal_distribution<double> distribution(0.0,sigma);

  arma::vec force = arma::zeros<arma::vec>(dimensions);
  //  force = -gamma*system[j].mass*system[j].velocity;
  for(int i=0;i<dimensions;i++){
    force(i) += distribution(generator);
  }
  return force;
p}
arma::vec particle_system::electrostatic_force(){
  int index=0;
  int index_r;
  int dim=this->dimensions+1;
  double force_ij=0.0;
  double r3=0.0;
  arma::vec relcoord = relative_coordinate();
  arma::vec force = arma::zeros<arma::vec>(particle_count*dimensions);
  for(int i=0;i<particle_count;i++){
    for(int j=0;j<i;j++){
      index_r = (i*(i-1)/2 + j)*dim + dim-1;
      printf("%d %f\n",index_r,relcoord(index_r));
      r3 = relcoord(index_r)*relcoord(index_r)*relcoord(index_r);
      for(int k=0;k<dimensions;k++){
	index = (i*(i-1)/2 + j)*(dim) + k;
	force_ij = system[i].charge*system[j].charge*relcoord(index)/r3;
	force(i*dimensions+k) += force_ij;
	force(j*dimensions+k) += -force_ij;
      }
    }
  }
  return force;      
}
void particle_system::verlet(int steps){
  double half_delta_t = delta_t/2.0;
  double half_delta_t_sq = delta_t*half_delta_t;
  double size = 10.0;
  
  
}
void particle_system::update(){
}
void particle_system::run(int steps){
}
void particle_system::print(){
  for(int i=0;i<particle_count;i++){
    printf("Particle%d:\n",i);
    this->system[i].print();
  }
}
