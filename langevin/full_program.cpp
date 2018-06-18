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
  particle(int Dimensions);
  
  void update(arma::vec position, arma::vec velocity);
};

class particle_system{
public:
  int particle_count; //number of particles
  int dimensions; //variable count
  int observables; //stats to be collected
  double size; //periodic boundary

  double delta_t; //time step + variance of langevin force
  double gamma; //drag factor
  double temp; //temperature
  double sigma; //sqrt(2*m*T*gamma*delta_t)

  std::vector<particle> system; //the system to be evolved
  arma::mat position; //position of system particles, used for verlet
  arma::mat velocity; //velocity of system particles, used for verlet
  arma::mat statistics; //average momentum, position

  particle_system(int Dimensions, double Size, double Delta_t, double Gamma, double Temp, int Particles); //auto-constructor  
  
  void add_periodic(particle new_particle); //add particle & convert positions to periodic boundary
  arma::vec relative_coordinate_periodic(); //calc relative coordinates to periodic boundary
  arma::vec langevin_force(); //add langevin force
  arma::vec electrostatic_force_periodic(); //basic repulsive force to periodic boundary

  void verlet_periodic(int steps, bool repulsion, int print_checks); //evolve one time step using periodic

  void statistics_run();
};




/******************** 
  CLASS FUNCTIONS
 ********************/
/******************** 
  particle
 ********************/
particle::particle(){
  this->dimensions = 0;
  this->mass = 0;
  this->charge = 0;
  this->position = arma::zeros<arma::vec>(0);
  this->velocity = arma::zeros<arma::vec>(0);
}
particle::particle(int Dimensions){
  this->dimensions = Dimensions;
  this->mass = 1.0;
  this->charge = 1.0;
  this->position = arma::zeros<arma::vec>(dimensions+1);
  this->velocity = arma::zeros<arma::vec>(dimensions+1);
}
void particle::update(arma::vec Position, arma::vec Velocity){
  this->position = Position;
  this->velocity = Velocity;
}

/******************** 
  particle_system
 ********************/
particle_system::particle_system(int Dimensions, double Size, double Delta_t, double Gamma, double Temp, int Particles){
  int particles=Particles;
  this->dimensions = Dimensions;
  this->size = Size;
  this->delta_t = Delta_t;
  this->gamma = Gamma;
  this->temp = Temp;
  this->particle_count=0;
  this->sigma=sqrt(2*temp*dimensions*gamma*delta_t);
  this->observables=1+2*(1+Dimensions);
  this->statistics = arma::zeros<arma::mat>(2,observables);

  double x=0.0,v=0.0,pi=acos(-1);
  particle part(dimensions);
  for(int i=0;i<particles;i++){
    x=0.0;
    v=0.0;
    for(int j=0;j<dimensions;j++){
      part.position(j) = i*pi + j*pi*pi*pi;
      x += part.position(j)*part.position(j);
      part.velocity(j) = (j+1)*pow(-1,i);
      v += part.velocity(j)*part.velocity(j);
    }
    part.position(dimensions) = sqrt(x);
    part.velocity(dimensions) = sqrt(v);
    add_periodic(part);
  }
}
void particle_system::add_periodic(particle new_particle){
  double x=0.0,v=0.0;
  this->system.push_back(new_particle);
  for(int i=0;i<dimensions;i++){
    system[particle_count].position(i) = std::fmod(new_particle.position(i),size);
    while(system[particle_count].position(i) < 0){
      system[particle_count].position(i) += size;
    }
    while(system[particle_count].position(i) > size){
      system[particle_count].position(i) -= size;
    }
    system[particle_count].velocity(i) = new_particle.velocity(i);

    x += system[particle_count].position(i)*system[particle_count].position(i);
    v += system[particle_count].velocity(i)*system[particle_count].velocity(i);
  }
  system[particle_count].position(dimensions) = sqrt(x);
  system[particle_count].velocity(dimensions) = sqrt(v);
  
  this->position.resize(dimensions+1,particle_count+1);
  this->position.col(particle_count) = system[particle_count].position;
  this->velocity.resize(dimensions+1,particle_count+1);
  this->velocity.col(particle_count) = system[particle_count].velocity;
  this->particle_count++;
}
arma::vec particle_system::relative_coordinate_periodic(){
  int index=0;
  int dim=this->dimensions+1;
  double r2=0.0,xi=0.0,xj=0.0,magnitude=0.0,period_mag=0.0;
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
	  magnitude = fabs(xi-xj);
	  period_mag = size - magnitude;
	  if(magnitude < period_mag){
	    relcoord(index) = position(k,i) - position(k,j);
	  }else{
	    if(xi > xj){
	      relcoord(index) = xi - (size + xj);
	    }else{
	      relcoord(index) = (size + xi) - xj;
	    }
	  }
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
  std::normal_distribution<double> distribution(0.0,1.0);
  arma::vec force = arma::zeros<arma::vec>(particle_count*dimensions);
  for(int i=0;i<particle_count;i++){
    for(int j=0;j<dimensions;j++){
      force(i*dimensions+j) = -gamma*system[i].mass*velocity(j,i);
      force(i*dimensions+j) += sigma*distribution(generator);
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
void particle_system::verlet_periodic(int steps,bool repulsion, int print_checks){
  double half_delta_t = delta_t/2.0;
  double half_delta_tsq = delta_t*half_delta_t;
  double acc_jk,x,v;
  int dim = this->dimensions+1;
  arma::vec force = arma::zeros<arma::vec>(particle_count*dimensions);
  arma::vec relcoord = arma::zeros<arma::vec>(particle_count*(particle_count-1)/2*dim);
  statistics = arma::zeros<arma::mat>(2,observables);
  if(repulsion==true){
    force = electrostatic_force_periodic();
    force += langevin_force();
  } else{
    force = langevin_force();
  }
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  double sqrt_mass;
  std::default_random_engine generator (seed);
  std::normal_distribution<double> distribution(7.3,1.5);
  
  for(int i=0;i<steps;i++){
    /*
    for(int j=0;j<particle_count;j++){
      x=0.0;
      for(int k=0;k<dimensions;k++){
	acc_jk = force(j*dimensions+k)/system[j].mass;
	position(k,j) = position(k,j) + delta_t*velocity(k,j) + half_delta_tsq*acc_jk;
	position(k,j) = std::fmod(position(k,j),size);
	while(position(k,j) > size){
	  position(k,j) -= size;
	}
	while(position(k,j) < 0){
	  position(k,j) += size;
	}
	x += position(k,j)*position(k,j);
	velocity(k,j) = velocity(k,j) + half_delta_t*acc_jk;
      }
      position(dimensions,j) = sqrt(x);
    }
    if(repulsion==true){
      force = electrostatic_force_periodic();
      force += langevin_force();
    } else{
      force = langevin_force();
    }
    for(int j=0;j<particle_count;j++){
      v=0.0;
      for(int k=0;k<dimensions;k++){
	acc_jk = force(j*dimensions+k)/system[j].mass;
	velocity(k,j) += half_delta_t*acc_jk;
	v += velocity(k,j)*velocity(k,j);
      }
      velocity(dimensions,j) = sqrt(v);
    }
    */
    for(int j=0;j<particle_count;j++){
      for(int k=0;k<dimensions+1;k++){
	position(k,j) = distribution(generator);
	velocity(k,j) = distribution(generator);
      }
    }
    statistics_run();
    if((i%print_checks)==0){
      printf("----------\nstep: %d\n",i);
    }
  }
  for(int i=0;i<observables;i++){
    statistics(0,i) /= (double) steps;
    statistics(1,i) /= (double) steps;
    statistics(1,i) = sqrt(statistics(1,i)-statistics(0,i)*statistics(0,i));
  }
}
void particle_system::statistics_run(){
  double momentum=0.0,mass=0.0;
  arma::vec avg_positions = arma::zeros<arma::vec>(dimensions+1);
  arma::vec avg_velocity = arma::zeros<arma::vec>(dimensions+1);
  for(int i=0;i<particle_count;i++){
    mass = system[i].mass;
    momentum += mass*mass*velocity(dimensions,i)*velocity(dimensions,i);
    for(int j=0;j<dimensions+1;j++){
      avg_positions(j) += position(j,i);
      avg_velocity(j) += velocity(j,i);
    }
  }
  momentum /= (double) particle_count;
  statistics(0,0) += momentum;
  statistics(1,0) += momentum*momentum;
  for(int i=0;i<dimensions+1;i++){
    avg_positions(i) /= (double) particle_count;
    avg_velocity(i) /= (double) particle_count;
    statistics(0,i+1) += avg_positions(i);
    statistics(1,i+1) += avg_positions(i)*avg_positions(i);
    statistics(0,i+2+dimensions) += avg_velocity(i);
    statistics(1,i+2+dimensions) += avg_velocity(i)*avg_velocity(i);
  }
}
void particle_system::print(){
  for(int i=0;i<particle_count;i++){
    printf("Particle%d:\n",i);
    this->system[i].print();
  }
}
