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
  int particle_count; //number of particles
  int dimensions; //variable count
  double size; //periodic boundary

  double delta_t; //time step + variance of langevin force
  double gamma; //drag factor
  double temp; //temperature
  double sigma; //sqrt(2*m*T*gamma*delta_t)

  std::vector<particle> system; //the system to be evolved
  arma::mat position; //position of system particles, used for verlet
  arma::mat velocity; //velocity of system particles, used for verlet
  //  arma::vec (*force)(particle part1, particle part2); //choice of evolution force

  particle_system(int Dimensions, double Size, double Delta_t, double Gamma, double Temp); //constructor
  
  void add(particle new_particle); //add particle
  void add_periodic(particle new_particle, bool True); //add particle & convert positions to periodic boundary
  arma::vec relative_coordinate(); //calc relative coordinates (i-j)
  arma::vec relative_coordinate_periodic(); //calc relative coordinates to periodic boundary
  arma::vec langevin_force(); //add langevin force
  arma::vec electrostatic_force(); //basic repulsive force
  arma::vec electrostatic_force_periodic(); //basic repulsive force to periodic boundary

  arma::mat verlet(int steps); //evolve one time step using verlet
  arma::mat verlet_periodic(int steps); //evolve one time step using periodic

  void update(); //update system info
  void run(int steps); //run system 'steps' times
  void print(); //print system information
};

class map_system{
public:
  int count;
  double size;
  arma::mat map;

  map_system(double Size);

  void add(double min, double max);
  double map_x_to_y(double x, int n);
};

/******************** 
  MAIN FUNCTION
 ********************/
std::ofstream ofile;

int main(int argc, char* argv[]){
  double size=10.0;
  double delta_t=0.05;
  double gamma=0.2;
  double temp=80;
  bool periodic=true;

  int dimensions=2;
  int particles=10;
  particle_system system(dimensions, size, delta_t, gamma, temp);
  arma::vec x = arma::zeros<arma::vec>(dimensions);
  arma::vec v = arma::zeros<arma::vec>(dimensions);
  arma::vec force_sum = arma::zeros<arma::vec>(dimensions);
  particle part1(dimensions,1,1,x,v);
  for(int i=0;i<particles;i++){
    part1.mass = 1.0;
    part1.charge = 1.0;
    for(int j=0;j<dimensions;j++){
      part1.position(j) = (i+j)*8.12312+7.234123;
      part1.velocity(j) = j*pow(-1,i);
    }
    system.add_periodic(part1,periodic);
  }

  system.position.print();
  system.velocity.print();
  arma::mat result = system.verlet_periodic(100);

  ofile.open("result.dat");
  for(int i=0;i<particles;i++){
    for(int j=0;j<dimensions;j++){
      ofile << std::setw(15) << result(j,i);
    }
    ofile << std::endl;
  }
  ofile.close();

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

particle_system::particle_system(int Dimensions, double Size, double Delta_t, double Gamma, double Temp){
  this->dimensions = Dimensions;
  this->size = Size;
  this->delta_t = Delta_t;
  this->gamma = Gamma;
  this->temp = Temp;
  this->particle_count=0;
  this->sigma=sqrt(2*Temp*gamma*delta_t);
}
void particle_system::add(particle new_particle){
  this->system.push_back(new_particle);
  this->particle_count++;
  this->position.resize(dimensions,particle_count);
  this->position.col(particle_count-1) = new_particle.position;
  this->velocity.resize(dimensions,particle_count);
  this->position.col(particle_count-1) = new_particle.velocity;
}
void particle_system::add_periodic(particle new_particle, bool True){
  this->system.push_back(new_particle);
  for(int i=0;i<dimensions;i++){
    system[particle_count].position(i) = std::fmod(new_particle.position(i),size);
    if(system[particle_count].position(i) < 0){
      system[particle_count].position(i) += size;
    }else if(system[particle_count].position(i) > size){
      system[particle_count].position(i) -= size;
    }
  }
  this->position.resize(dimensions,particle_count+1);
  this->position.col(particle_count) = system[particle_count].position;
  this->velocity.resize(dimensions,particle_count+1);
  this->velocity.col(particle_count) = new_particle.velocity;
  this->particle_count++;
}
arma::vec particle_system::relative_coordinate_periodic(){
  int index=0;
  int dim=this->dimensions+1;
  double r2=0.0,pi=0.0,pj=0.0,magnitude=0.0,period_mag=0.0;
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
	  pi = position(k,i);
	  pj = position(k,j);
	  magnitude = fabs(pi-pj);
	  period_mag = size - magnitude;
	  if(magnitude < period_mag){
	    relcoord(index) = position(k,i) - position(k,j);
	  }else{
	    if(pi > pj){
	      relcoord(index) = pi - (size + pj);
	    }else{
	      relcoord(index) = (size + pi) - pj;
	    }
	  }
	  r2 += relcoord(index)*relcoord(index);
	}
      }
    }       
  }
  return relcoord;
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
  double mass;
  std::default_random_engine generator (seed);
  std::normal_distribution<double> distribution(0.0,sqrt(delta_t));
  arma::vec force = arma::zeros<arma::vec>(particle_count*dimensions);
  for(int i=0;i<particle_count;i++){
    mass = sqrt(system[i].mass);
    for(int j=0;j<dimensions;j++){
      force(i*dimensions+j) = -gamma*system[i].mass*system[i].velocity(j);
      force(i*dimensions+j) += sigma*mass*distribution(generator);
    }
  }
  return force;
}
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
arma::vec particle_system::electrostatic_force_periodic(){
  int index=0;
  int index_r;
  int dim=this->dimensions+1;
  double force_ij=0.0;
  double r3=0.0;
  arma::vec relcoord = relative_coordinate_periodic();
  arma::vec force = arma::zeros<arma::vec>(particle_count*dimensions);
  for(int i=0;i<particle_count;i++){
    for(int j=0;j<i;j++){
      index_r = (i*(i-1)/2 + j)*dim + dim-1;
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
arma::mat particle_system::verlet(int steps){
  double half_delta_t = delta_t/2.0;
  double half_delta_tsq = delta_t*half_delta_t;
  double acc_jk;
  int dim = this->dimensions+1;
  arma::vec force = arma::zeros<arma::vec>(particle_count*dimensions);
  arma::vec relcoord = arma::zeros<arma::vec>(particle_count*(particle_count-1)/2*dim);
  force = electrostatic_force();
  force += langevin_force();

  for(int i=0;i<steps;i++){
    for(int j=0;j<particle_count;j++){
      for(int k=0;k<dimensions;k++){
	acc_jk = force(j*dimensions+k)/system[j].mass;
	position(k,j) = position(k,j) + delta_t*velocity(k,j) + half_delta_tsq*acc_jk;
	velocity(k,j) = velocity(k,j) + half_delta_t*acc_jk;
      }
    }
    force = electrostatic_force();
    force += langevin_force();
    for(int j=0;j<particle_count;j++){
      for(int k=0;k<dimensions;k++){
	acc_jk = force(j*dimensions+k)/system[j].mass;
	velocity(k,j) += half_delta_t*acc_jk;
      }
    }
    printf("step %d:\n",i);
    position.print();
    velocity.print();
    printf("-------------------\n");
  }
  return position;
}
arma::mat particle_system::verlet_periodic(int steps){
  double half_delta_t = delta_t/2.0;
  double half_delta_tsq = delta_t*half_delta_t;
  double acc_jk;
  int dim = this->dimensions+1;
  arma::vec force = arma::zeros<arma::vec>(particle_count*dimensions);
  arma::vec relcoord = arma::zeros<arma::vec>(particle_count*(particle_count-1)/2*dim);
  force = electrostatic_force_periodic();
  force += langevin_force();

  for(int i=0;i<steps;i++){
    for(int j=0;j<particle_count;j++){
      for(int k=0;k<dimensions;k++){
	acc_jk = force(j*dimensions+k)/system[j].mass;
	position(k,j) = position(k,j) + delta_t*velocity(k,j) + half_delta_tsq*acc_jk;
	position(k,j) = std::fmod(position(k,j),size);
	if(position(k,j) > size){
	  position(k,j) -= size;
	} else if(position(k,j) < 0){
	  position(k,j) += size;
	}
	velocity(k,j) = velocity(k,j) + half_delta_t*acc_jk;
      }
    }
    force = electrostatic_force_periodic();
    force += langevin_force();
    for(int j=0;j<particle_count;j++){
      for(int k=0;k<dimensions;k++){
	acc_jk = force(j*dimensions+k)/system[j].mass;
	velocity(k,j) += half_delta_t*acc_jk;
      }
    }
  printf("step %d:\n",i);
  position.print();
  velocity.print();
  printf("-------------------\n");
  }
  return position;
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

/******************** 
  map
 ********************/
map_system::map_system(double Size){
  this->size = Size;  
}
void map_system::add(double min, double max){
  map.resize(3,count);
  map(0,count)=min;
  map(1,count)=max;
  map(2,count)=(max-min)/size;
  count++;
}
double map_system::map_x_to_y(double x, int n){
  return x*map(2,n) + map(0,n);
}
