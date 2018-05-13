#ifndef PROJECT3_LIBRARY_H
#define PROJECT3_LIBRARY_H

#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>

using namespace std;

//massivebody contains the mass, positions, and velocities
//a constructor with a "primary" massivebody is given to create
//	an orbiting mass around said "primary"
class massivebody{
	public:
		double mass;
		double position[4];
		double velocity[4];
		massivebody();
		massivebody(double, double, double, double, double, double, double);
		massivebody(massivebody, double, double, double, double, double, double, double);
};

//massivesystem contains all necessary components to solve for evolution
//massivebody_count tracks the size of massivebody* system
//the total mass and Center of Mass positions & velocities are given via
//	mass_total, composition[4], comvelocity[4]
//after defining a massive system, 'add' resizes massivebody* and tacks the
//	new massivebody onto the end
class massivesystem{
	public:
		int massivebody_count;
		massivebody* system;
		double mass_total;
		double composition[4];
		double comvelocity[4];

		massivesystem();
		~massivesystem();

		void add(massivebody guy);
		void update(double**&, int);
		void relpositions(double**&);
		void forces(double*&, const double**&);
		void initialize(double**&, double*&, int);
		void verlet(double**&, int, double, int);
		void RK4(double**&, int, double, int);
		void print();
};
class vector{
	public:
		int length;
		double* array;
};

inline void array_alloc(double*& a, int length);
inline void array_delete(double*& a);
inline void array_resize(double*&, int, int);
inline void matrix_alloc(double**&, int, int);
inline void matrix_delete(double**&, int);
inline void matrix_resize(double**&, int, int, int, int);
inline void threeDarray_alloc(double***& a, int d1, int d2, int d3);
inline void threeDarray_delete(double***& a, int d1, int d2);
inline void rel_position(double**&, double**&, int, int);
inline void rel_coord(double**&, double**&, int, int);
inline void rel_angmomentum(double**&, double**&, int);
inline void gravityforces(double**&, double**&, double*&, int, int);
inline void gravityforces_relcorrection(double**&, double**&, double**&, double*&, int, int);
void verlet(double**&, double*&, int, int, double, int);
void verlet_relcor(double**&, double*&, int, int, double, int);
void helionstates(double**&, double**&, double**&, int, int);
void helionstates_dynamic(double**& aphelion, int& aphguess, double**& perihelion, int& periguess, double**& output, int steps, int body, int range, int dec_or_inc);

/**************************
Begin function definitions
***************************/


/*********************
massivebody functions
*********************/
massivebody::massivebody(){
	mass = 0.0;
	position[0] = 0.0;
	position[1] = 0.0;
	position[2] = 0.0;
	position[3] = 0.0;
	velocity[0] = 0.0;
	velocity[1] = 0.0;
	velocity[2] = 0.0;
	velocity[3] = 0.0;
}
massivebody::massivebody(double Mass, double x, double y, double z, double vx, double vy, double vz){
	mass = Mass;
	position[0] = x;
	position[1] = y;
	position[2] = z;
	position[3] = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
	velocity[0] = vx;
	velocity[1] = vy;
	velocity[2] = vz;
	velocity[3] = sqrt(vx*vx + vy*vy + vz*vz);
}
massivebody::massivebody(massivebody primary, double Mass, double x, double y, double z, double vx, double vy, double vz){
	mass = Mass;
	position[0] = primary.position[0] + x;
	position[1] = primary.position[1] + y;
	position[2] = primary.position[2] + z;
	position[3] = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
	velocity[0] = primary.velocity[0] + vx;
	velocity[1] = primary.velocity[1] + vy;
	velocity[2] = primary.velocity[2] + vz;
	velocity[3] = sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
}

/***********************
massivesystem functions
***********************/
massivesystem::massivesystem(){
	massivebody_count = 0;
	system = nullptr;
	mass_total = 0;
	for(int i=0;i<4;i++){
		composition[i]=0.0;
		comvelocity[i]=0.0;
	}
}
massivesystem::~massivesystem(){
	delete[] system;
}
void massivesystem::add(massivebody guy){
	double mass_old = mass_total;
	massivebody* temp = new massivebody[massivebody_count+1];
	for(int i=0;i<massivebody_count;i++){
		temp[i].mass=system[i].mass;
		temp[i].position[0]=system[i].position[0];
		temp[i].position[1]=system[i].position[1];
		temp[i].position[2]=system[i].position[2];
		temp[i].position[3]=system[i].position[3];
		temp[i].velocity[0]=system[i].velocity[0];
		temp[i].velocity[1]=system[i].velocity[1];
		temp[i].velocity[2]=system[i].velocity[2];
		temp[i].velocity[3]=system[i].velocity[3];
	}
	temp[massivebody_count].mass=guy.mass;
	temp[massivebody_count].position[0]=guy.position[0];
	temp[massivebody_count].position[1]=guy.position[1];
	temp[massivebody_count].position[2]=guy.position[2];
	temp[massivebody_count].position[3]=guy.position[3];
	temp[massivebody_count].velocity[0]=guy.velocity[0];
	temp[massivebody_count].velocity[1]=guy.velocity[1];
	temp[massivebody_count].velocity[2]=guy.velocity[2];
	temp[massivebody_count].velocity[3]=guy.velocity[3];

	massivebody_count++;

	delete[] system;
	mass_total += guy.mass;
	for(int i=0;i<3;i++){
		composition[i] = (composition[i]*mass_old + guy.position[i]*guy.mass)/mass_total;
		comvelocity[i] = (comvelocity[i]*mass_old + guy.velocity[i]*guy.mass)/mass_total;
	}
	composition[3] = sqrt(composition[0]*composition[0] + composition[1]*composition[1] + composition[2]*composition[2]);
	comvelocity[3] = sqrt(comvelocity[0]*comvelocity[0] + comvelocity[1]*comvelocity[1] + comvelocity[2]*comvelocity[2]);
	system = temp;
}
void massivesystem::update(double**& output, int i){
	int j;
	mass_total = 0.0;
	for(j=0;j<4;j++){
		composition[j] = 0.0;
		comvelocity[j] = 0.0;
	}

	for(j=0;j<massivebody_count;j++){
		system[j].position[0] = output[i][8*j+1];
		system[j].position[1] = output[i][8*j+2];
		system[j].position[2] = output[i][8*j+3];
		system[j].position[3] = output[i][8*j+4];
		system[j].velocity[0] = output[i][8*j+5];
		system[j].velocity[1] = output[i][8*j+6];
		system[j].velocity[2] = output[i][8*j+7];
		system[j].velocity[3] = output[i][8*j+8];
		mass_total += system[j].mass;
	}
	for(j=0;j<massivebody_count;j++){
		composition[0] += system[j].position[0]*system[j].mass;
		composition[1] += system[j].position[1]*system[j].mass;
		composition[2] += system[j].position[2]*system[j].mass;
		comvelocity[0] += system[j].velocity[0]*system[j].mass;
		comvelocity[1] += system[j].velocity[1]*system[j].mass;
		comvelocity[2] += system[j].velocity[2]*system[j].mass;
	}
	comvelocity[0] /= mass_total;
	comvelocity[1] /= mass_total;
	comvelocity[2] /= mass_total;

	composition[3] = sqrt(composition[0]*composition[0] + composition[1]*composition[1] + composition[2]*composition[2]);
	comvelocity[3] = sqrt(comvelocity[0]*comvelocity[0] + comvelocity[1]*comvelocity[1] + comvelocity[2]*comvelocity[2]);
}
void massivesystem::print(){
	cout << "mass=" << setw(15) << mass_total << endl;
	cout << "body_count=" << setw(15) << massivebody_count << endl;
	for(int i=0;i<4;i++){
		cout << "pos[" << i<< "]=" << setw(15) << composition[i] << endl;
	}
	for(int i=0;i<4;i++){
		cout << "vel[" << i<< "]=" << setw(15) << comvelocity[i] << endl;
	}
}
void massivesystem::relpositions(double**& relcoord){
	int i, j;	
	int bodies = massivebody_count;
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies;j++){
			if(i!=j){
				relcoord[i][j*4] = system[i].position[0] - system[j].position[0];
				relcoord[i][j*4+1] = system[i].position[1] - system[j].position[1];
				relcoord[i][j*4+2] = system[i].position[2] - system[j].position[2];
				relcoord[i][j*4+3] = sqrt(relcoord[i][j*4]*relcoord[i][j*4] + relcoord[i][j*4+1]*relcoord[i][j*4+1]+relcoord[i][j*4+2]*relcoord[i][j*4+2]);
			}
		}
	}
}
void massivesystem::forces(double*& force, const double**& relcoord){
	int j, k;
	int bodies = massivebody_count;
	double rcubed, fourpisq;
	fourpisq = 4.0*M_PI*M_PI;
	for(j=0;j<bodies;j++){
		force[j*3] = 0.0;
		force[j*3+1] = 0.0;
		force[j*3+2] = 0.0;
	}
	for(j=0;j<bodies;j++){
		for(k=0;k<bodies;k++){
			if(j!=k){
				rcubed = pow(relcoord[j][k*4+3],3.0);
				force[j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
				force[j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
				force[j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
			}
		}
		force[j*3] *= fourpisq;
		force[j*3+1] *= fourpisq;
		force[j*3+2] *= fourpisq;
	}
}
void massivesystem::initialize(double**& output, double*& mass, int steps){
	int i, j;
	int bodies = massivebody_count;
	output = new double*[steps+1];
	mass = new double[bodies];
	for(int i=0;i<steps+1;i++){
		output[i] = new double[bodies*8+1];
	}
	for(i=1;i<steps+1;i++){
		for(j=0;j<bodies;j++){
			output[i][j*8+1] = 0;
			output[i][j*8+2] = 0;
			output[i][j*8+3] = 0;
			output[i][j*8+4] = 0;
			output[i][j*8+5] = 0;
			output[i][j*8+6] = 0;
			output[i][j*8+7] = 0;
			output[i][j*8+8] = 0;
		}
	}
	for(j=0;j<bodies;j++){
		mass[j] = system[j].mass;
		output[0][j*8+1] = system[j].position[0];
		output[0][j*8+2] = system[j].position[1];
		output[0][j*8+3] = system[j].position[2];
		output[0][j*8+4] = system[j].position[3];
		output[0][j*8+5] = system[j].velocity[0];
		output[0][j*8+6] = system[j].velocity[1];
		output[0][j*8+7] = system[j].velocity[2];
		output[0][j*8+8] = system[j].velocity[3];
	}
}


/***************************
class function solvers,
RK4 and Verlet
***************************/
void massivesystem::RK4(double**& output, int steps, double tmax, int sun_choice){
	int i, j, k, m, n, bodies;
	bodies = massivebody_count;
	double h, halfh, halfhsq, rcubed, fourpisq;
	double* k2, *k3;
	h = (tmax-0.0)/((double) (steps));
	halfh = 0.5*h;
	halfhsq = halfh*h;	
	fourpisq = 4.0*M_PI*M_PI;

	//placeholders arrays
	k2 = new double[bodies*6];
	k3 = new double[bodies*6];
	for(i=0;i<bodies*6;i++){
		k2[i] = 0;
		k3[i] = 0;
	}
	//relcoord to provide relative coordinates, x, y, z, r, of body i with respect to body j
	double** relcoord = new double*[bodies];
	for(i=0;i<bodies;i++){
		relcoord[i] = new double[bodies*4];
	}
	//force provides fx,fy,fz for each body at step i and i+1
	double** force = new double*[2];
	for(i=0;i<2;i++){
		force[i] = new double[bodies*3];
	}
	//initialize all matrices to zero
	for(i=0;i<steps+1;i++){
		for(j=0;j<bodies*8+1;j++){
			output[i][j]=0.0;
		}
	}
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies*4;j++){
			relcoord[i][j]=0.0;
		}
	}
	for(i=0;i<2;i++){
		for(j=0;j<bodies*3;j++){
			force[i][j]=0.0;
		}
	}

	//reinitialize relcoord to initial conditions
	//technically double calculating--consider revising
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies;j++){
			if(i!=j){
				relcoord[i][j*4] = system[i].position[0] - system[j].position[0];
				relcoord[i][j*4+1] = system[i].position[1] - system[j].position[1];
				relcoord[i][j*4+2] = system[i].position[2] - system[j].position[2];
				relcoord[i][j*4+3] = sqrt(relcoord[i][j*4]*relcoord[i][j*4] + relcoord[i][j*4+1]*relcoord[i][j*4+1]+relcoord[i][j*4+2]*relcoord[i][j*4+2]);
			}
		}
	}
	//initialize forces (k1v)
	for(j=0;j<bodies;j++){
		for(k=0;k<bodies;k++){
			if(j!=k){
				rcubed = pow(relcoord[j][k*4+3],3.0);
				force[0][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
				force[0][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
				force[0][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
			}
		}
		force[0][j*3] *= fourpisq;
		force[0][j*3+1] *= fourpisq;
		force[0][j*3+2] *= fourpisq;
	}
	//initialize initial conditions for output
	for(j=0;j<bodies;j++){
		output[0][j*8+1] = system[j].position[0];
		output[0][j*8+2] = system[j].position[1];
		output[0][j*8+3] = system[j].position[2];
		output[0][j*8+4] = system[j].position[3];
		output[0][j*8+5] = system[j].velocity[0];
		output[0][j*8+6] = system[j].velocity[1];
		output[0][j*8+7] = system[j].velocity[2];
		output[0][j*8+8] = system[j].velocity[3];
	}

	//Begin RK4 loop
	for(i=0;i<steps;i++){
		output[i+1][0] = output[i][0] + h;
		//calculate first positions [1+1/2]
		// (using k1x)
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+1] = output[i][j*8+1] + halfh*output[i][j*8+5];
			output[i+1][j*8+2] = output[i][j*8+2] + halfh*output[i][j*8+6];
			output[i+1][j*8+3] = output[i][j*8+3] + halfh*output[i][j*8+7];
		}
		//new relative positions at [i+1/2]
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
					relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
					relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//new forces with relative positions
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
			//k2v
			k2[j*6+3] = force[1][j*3];
			k2[j*6+4] = force[1][j*3+1];
			k2[j*6+5] = force[1][j*3+2];
		}
		//k2x using i+1/2 positions
		for(j=0;j<bodies;j++){
			k2[j*6] = output[i][j*8+5] + halfh*force[1][j*3];
			k2[j*6+1] = output[i][j*8+6] + halfh*force[1][j*3+1];
			k2[j*6+2] = output[i][j*8+7] + halfh*force[1][j*3+2];
		}
		//new temp positions at i+1/2 with k2x
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+1] = output[i][j*8+1] + halfh*k2[j*6];
			output[i+1][j*8+2] = output[i][j*8+2] + halfh*k2[j*6+1];
			output[i+1][j*8+3] = output[i][j*8+3] + halfh*k2[j*6+2];
		}
		//new relative positions at [i+1/2] with k2x
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
					relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
					relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//new forces with relative positions
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
			//k3v
			k3[j*6+3] = force[1][j*3];
			k3[j*6+4] = force[1][j*3+1];
			k3[j*6+5] = force[1][j*3+2];
		}
		//k3x using i+1/2 positions
		for(j=0;j<bodies;j++){
			k3[j*6] = output[i][j*8+5] + halfh*force[1][j*3];
			k3[j*6+1] = output[i][j*8+6] + halfh*force[1][j*3+1];
			k3[j*6+2] = output[i][j*8+7] + halfh*force[1][j*3+2];
		}
		//new temp positions at i+1 using k3x
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+1] = output[i][j*8+1] + h*k3[j*6];
			output[i+1][j*8+2] = output[i][j*8+2] + h*k3[j*6+1];
			output[i+1][j*8+3] = output[i][j*8+3] + h*k3[j*6+2];
		}
		//new relative positions at [i+1] with k3x
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
					relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
					relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//new forces with relative positions
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
		}
		//k4x using i+1/2 positions
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+5] = output[i][j*8+5] + h*force[1][j*3];
			output[i+1][j*8+6] = output[i][j*8+6] + h*force[1][j*3+1];
			output[i+1][j*8+7] = output[i][j*8+7] + h*force[1][j*3+2];
		}
		//final positions at i+1
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+1] = output[i][j*8+1] + (h/6.0)*(output[i][j*8+5] + 2.0*k2[j*6] + 2.0*k3[j*6] + output[i+1][j*8+5]);
			output[i+1][j*8+2] = output[i][j*8+2] + (h/6.0)*(output[i][j*8+6] + 2.0*k2[j*6+1] + 2.0*k3[j*6+1] + output[i+1][j*8+6]);
			output[i+1][j*8+3] = output[i][j*8+3] + (h/6.0)*(output[i][j*8+7] + 2.0*k2[j*6+2] + 2.0*k3[j*6+2] + output[i+1][j*8+7]);
			output[i+1][j*8+4] = sqrt(output[i+1][j*8+1]*output[i+1][j*8+1] + output[i+1][j*8+2]*output[i+1][j*8+2] + output[i+1][j*8+3]*output[i+1][j*8+3]);
		}	
		//final velocities at i+1
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+5] = output[i][j*8+5] + (h/6.0)*(force[0][j*3] + 2.0*k2[j*6+3] + 2.0*k3[j*6+3] + force[1][j*3]);
			output[i+1][j*8+6] = output[i][j*8+6] + (h/6.0)*(force[0][j*3+1] + 2.0*k2[j*6+4] + 2.0*k3[j*6+4] + force[1][j*3+1]);
			output[i+1][j*8+7] = output[i][j*8+7] + (h/6.0)*(force[0][j*3+2] + 2.0*k2[j*6+5] + 2.0*k3[j*6+5] + force[1][j*3+2]);
			output[i+1][j*8+8] = sqrt(output[i+1][j*8+5]*output[i+1][j*8+5] + output[i+1][j*8+6]*output[i+1][j*8+6] + output[i+1][j*8+7]*output[i+1][j*8+7]);
		}
		//final relative positions at i+1
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
					relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
					relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//final forces with relative positions
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
		}
		for(j=0;j<bodies;j++){
			force[0][j*3] = force[1][j*3];
			force[0][j*3+1] = force[1][j*3+1];
			force[0][j*3+2] = force[1][j*3+2];
		}

	}

	//delete dross--keep output
	for(i=0;i<bodies;i++){
		delete[] relcoord[i];
	}
	delete[] relcoord;
	for(i=0;i<2;i++){
		delete[] force[i];
	}
	delete[] force;
	delete[] k2;
	delete[] k3;
}
void massivesystem::verlet(double**& output, int steps, double tmax, int sun_choice){
	int i, j, k, m, n, bodies;
	bodies = massivebody_count;
	double h, halfh, halfhsq, rcubed, fourpisq;
	h = (tmax-0.0)/((double) (steps));
	halfh = 0.5*h;
	halfhsq = halfh*h;	
	fourpisq = 4.0*M_PI*M_PI;
	//output to provide t and x, y, z, vx, vy, vz for each body at each step i
	output = new double*[steps+1];	
	for(int i=0;i<steps+1;i++){
		output[i] = new double[bodies*8+1];
	}
	//relcoord to provide relative coordinates, x, y, z, r, of body i with respect to body j
	//technically double counting--consider revising
	double** relcoord = new double*[bodies];
	for(i=0;i<bodies;i++){
		relcoord[i] = new double[bodies*4];
	}
	//force provides fx,fy,fz for each body at step i and i+1
	double** force = new double*[2];
	for(i=0;i<2;i++){
		force[i] = new double[bodies*3];
	}
	//initialize all matrices to zero
	for(i=0;i<steps+1;i++){
		for(j=0;j<bodies*8+1;j++){
			output[i][j]=0.0;
		}
	}
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies*4;j++){
			relcoord[i][j]=0.0;
		}
	}
	for(i=0;i<2;i++){
		for(j=0;j<bodies*3;j++){
			force[i][j]=0.0;
		}
	}

	//reinitialize relcoord to initial conditions
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies;j++){
			if(i!=j){
				relcoord[i][j*4] = system[i].position[0] - system[j].position[0];
				relcoord[i][j*4+1] = system[i].position[1] - system[j].position[1];
				relcoord[i][j*4+2] = system[i].position[2] - system[j].position[2];
				relcoord[i][j*4+3] = sqrt(relcoord[i][j*4]*relcoord[i][j*4] + relcoord[i][j*4+1]*relcoord[i][j*4+1]+relcoord[i][j*4+2]*relcoord[i][j*4+2]);
			}
		}
	}
	//initialize forces
	for(j=0;j<bodies;j++){
		for(k=0;k<bodies;k++){
			if(j!=k){
				rcubed = pow(relcoord[j][k*4+3],3.0);
				force[0][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
				force[0][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
				force[0][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
			}
		}
		force[0][j*3] *= fourpisq;
		force[0][j*3+1] *= fourpisq;
		force[0][j*3+2] *= fourpisq;
	}
	//initialize initial conditions for output
	for(j=0;j<bodies;j++){
		output[0][j*8+1] = system[j].position[0];
		output[0][j*8+2] = system[j].position[1];
		output[0][j*8+3] = system[j].position[2];
		output[0][j*8+4] = system[j].position[3];
		output[0][j*8+5] = system[j].velocity[0];
		output[0][j*8+6] = system[j].velocity[1];
		output[0][j*8+7] = system[j].velocity[2];
		output[0][j*8+8] = system[j].velocity[3];
	}

	//begin solution for-loop
	for(i=0;i<steps;i++){
		//Set next time value
		output[i+1][0] = output[i][0] + h;
		//calculate x,y,z		
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+1] = output[i][j*8+1] + h*output[i][j*8+5] + halfhsq*force[0][j*3];
			output[i+1][j*8+2] = output[i][j*8+2] + h*output[i][j*8+6] + halfhsq*force[0][j*3+1];
			output[i+1][j*8+3] = output[i][j*8+3] + h*output[i][j*8+7] + halfhsq*force[0][j*3+2];
			output[i+1][j*8+4] = sqrt(output[i][j*8+1]*output[i][j*8+1] + output[i][j*8+2]*output[i][j*8+2] + output[i][j*8+3]*output[i][j*8+3]);
		}
		//new relcoord
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
					relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
					relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//calculate force_i+1
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
		}
		//calculate vx, vy, vz
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+5] = output[i][j*8+5] + halfh*(force[1][j*3] + force[0][j*3]);
			output[i+1][j*8+6] = output[i][j*8+6] + halfh*(force[1][j*3+1] + force[0][j*3+1]);
			output[i+1][j*8+7] = output[i][j*8+7] + halfh*(force[1][j*3+2] + force[0][j*3+2]);
			output[i+1][j*8+8] = sqrt(output[i+1][j*8+5]*output[i+1][j*8+5] + output[i+1][j*8+6]*output[i+1][j*8+6] + output[i+1][j*8+7]*output[i+1][j*8+7]);
		}
		//force_i for next loop is force_i+1 for this loop. No need for double computation.
		for(j=0;j<bodies;j++){
			force[0][j*3] = force[1][j*3];
			force[0][j*3+1] = force[1][j*3+1];
			force[0][j*3+2] = force[1][j*3+2];
		}
	}
	//delete dross--keep output
	for(i=0;i<bodies;i++){
		delete[] relcoord[i];
	}
	delete[] relcoord;
	for(i=0;i<2;i++){
		delete[] force[i];
	}
	delete[] force;

}


/******************************************************
	Begin function non-class functions
******************************************************/

/***************
Array functions
***************/
inline void array_alloc(double*& a, int length){
	int i;	
	a = new double[length];
	for(i=0;i<length;i++){
		a[i] = 0.0;
	}
}
inline void array_delete(double*& a){
	delete[] a;
}
inline void array_resize(double*& array, int oldsize, int newsize){
	double* temp = new double[newsize];
	for(int i=0;i<oldsize;i++){
		temp[i]=array[i];
	}

	delete[] array;
	array = temp;
}


/***************
Matrix functions
***************/
inline void matrix_alloc(double**& a, int rows, int columns){
	int i, j;	
	a = new double*[rows];
	for(i=0;i<rows;i++){
		a[i] = new double[columns];
	}
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			a[i][j] = 0.0;
		}
	}
}
inline void matrix_delete(double**& a, int rows){
	for(int i=0;i<rows;i++){
		delete a[i];
	}
	delete[] a;
}
inline void matrix_resize(double**& matrix, int oldrows, int oldcol, int newrows, int newcol){
	int i, j;	
	double** temp = new double*[newrows];
	for(i=0;i<newrows;i++){
		temp[i] = new double[newcol];
	}
	for(i=0;i<oldrows;i++){
		for(j=0;j<oldcol;j++){
			temp[i][j]=matrix[i][j];
		}
	}

	for(i=0;i<oldrows;i++){
		delete[] matrix[i];
	}
	delete[] matrix;
	matrix = temp;
}


/*****************
3D-Array functions
*****************/
inline void threeDarray_alloc(double***& a, int d1, int d2, int d3){
	int i, j, k;	
	a = new double**[d1];
	for(i=0;i<d1;i++){
		a[i] = new double*[d2];
	}
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			a[i][j] = new double[d3];
		}
	}
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				a[i][j][k] = 0.0;
			}
		}
	}
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				a[i][j][k] = 0.0;
			}
		}
	}
}
inline void threeDarray_delete(double***& a, int d1, int d2){
	int i, j;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			delete[] a[i][j];
		}
	}
	for(i=0;i<d1;i++){
		delete[] a[i];
	}
	delete[] a;
}


/***************************
force calculation functions
***************************/
inline void rel_position(double**& relcoord, double**& output, int i, int bodies){
	int m, n;
	for(m=0;m<bodies;m++){
		for(n=0;n<bodies;n++){
			if(m!=n){
				relcoord[m][n*4] = output[i][m*8+1] - output[i][n*8+1];
				relcoord[m][n*4+1] = output[i][m*8+2] - output[i][n*8+2];
				relcoord[m][n*4+2] = output[i][m*8+3] - output[i][n*8+3];
				relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
			}
		}
	}	
}
inline void rel_coord(double**& relcoord, double**& output, int i, int bodies){
	int m, n;
	for(m=0;m<bodies;m++){
		for(n=0;n<bodies;n++){
			if(m!=n){
				relcoord[m][n*8] = output[i][m*8+1] - output[i][n*8+1];
				relcoord[m][n*8+1] = output[i][m*8+2] - output[i][n*8+2];
				relcoord[m][n*8+2] = output[i][m*8+3] - output[i][n*8+3];
				relcoord[m][n*8+3] = sqrt(relcoord[m][n*8]*relcoord[m][n*8] + relcoord[m][n*8+1]*relcoord[m][n*8+1]+relcoord[m][n*8+2]*relcoord[m][n*8+2]);
				relcoord[m][n*8+4] = output[i][m*8+5] - output[i][n*8+5];
				relcoord[m][n*8+5] = output[i][m*8+6] - output[i][n*8+6];
				relcoord[m][n*8+6] = output[i][m*8+7] - output[i][n*8+7];
				relcoord[m][n*8+7] = sqrt(relcoord[m][n*8+4]*relcoord[m][n*8+4] + relcoord[m][n*8+5]*relcoord[m][n*8+5]+relcoord[m][n*8+6]*relcoord[m][n*8+6]);
			}
		}
	}	
}
inline void rel_angmomentum(double**& relangmomentum, double**& relcoord, int bodies){
	int m, n;
	for(m=0;m<bodies;m++){
		for(n=0;n<bodies;n++){
			if(m!=n){
				relangmomentum[m][n*4] = relcoord[m][n*8+1]*relcoord[m][n*8+6] - relcoord[m][n*8+2]*relcoord[m][n*8+5];
				relangmomentum[m][n*4+1] = relcoord[m][n*8+2]*relcoord[m][n*8+4] - relcoord[m][n*8]*relcoord[m][n*8+6];
				relangmomentum[m][n*4+2] = relcoord[m][n*8]*relcoord[m][n*8+5] - relcoord[m][n*8+1]*relcoord[m][n*8+4];
				relangmomentum[m][n*4+3] = sqrt(relangmomentum[m][n*4]*relangmomentum[m][n*4] + relangmomentum[m][n*4+1]*relangmomentum[m][n*4+1] + relangmomentum[m][n*4+2]*relangmomentum[m][n*4+2]);
			}
		}
	}
}
inline void gravityforces(double**& force, double**& relcoord, double*& mass, int i, int bodies){
	int j, k;
	double rcubed, fourpisq;
	fourpisq = 4.0*M_PI*M_PI;
	for(j=0;j<bodies;j++){
		force[i][j*3] = 0.0;
		force[i][j*3+1] = 0.0;
		force[i][j*3+2] = 0.0;
	}
	for(j=0;j<bodies;j++){
		for(k=0;k<bodies;k++){
			if(j!=k){
				rcubed = pow(relcoord[j][k*4+3],3.0);
				force[i][j*3] += -mass[k]*relcoord[j][k*4]/rcubed;
				force[i][j*3+1] += -mass[k]*relcoord[j][k*4+1]/rcubed;
				force[i][j*3+2] += -mass[k]*relcoord[j][k*4+2]/rcubed;
			}
		}
		force[i][j*3] *= fourpisq;
		force[i][j*3+1] *= fourpisq;
		force[i][j*3+2] *= fourpisq;
	}
}
inline void gravityforces_relcorrection(double**& force, double**& relangmomentum, double**& relcoord, double*& mass, int i, int bodies){
	int j, k;
	double rsq, rcubed, fourpisq, relcor, threeovercsq, fullterm;
	threeovercsq = 3.0/(63197.8*63197.8);
	fourpisq = 4.0*M_PI*M_PI;
	for(j=0;j<bodies;j++){
		force[i][j*3] = 0.0;
		force[i][j*3+1] = 0.0;
		force[i][j*3+2] = 0.0;
	}
	for(j=0;j<bodies;j++){
		for(k=0;k<bodies;k++){
			if(j!=k){
				rsq = pow(relcoord[j][k*8+3],2.0);
				relcor = (threeovercsq/rsq)*(relangmomentum[j][k*4+3]*relangmomentum[j][k*4+3]);
				fullterm = (1.0 + relcor)/(rsq*relcoord[j][k*8+3]);
				force[i][j*3] += -mass[k]*relcoord[j][k*8]*fullterm;
				force[i][j*3+1] += -mass[k]*relcoord[j][k*8+1]*fullterm;
				force[i][j*3+2] += -mass[k]*relcoord[j][k*8+2]*fullterm;
			}
		}
		force[i][j*3] *= fourpisq;
		force[i][j*3+1] *= fourpisq;
		force[i][j*3+2] *= fourpisq;
	}
}


/***************************
Non-class function Verlet
and relativistic Verlet
***************************/
void verlet(double**& output, double*& mass, int bodies, int steps, double tmax, int sun_choice){
	int i, j, k, m, n;
	double h, halfh, halfhsq, rcubed, fourpisq;
	h = (tmax-0.0)/((double) (steps));
	halfh = 0.5*h;
	halfhsq = halfh*h;	
	fourpisq = 4.0*M_PI*M_PI;

	//relcoord to provide relative coordinates, x, y, z, r, of body i with respect to body j
	//technically double calculating--consider revising
	double** relcoord = new double*[bodies];
	for(i=0;i<bodies;i++){
		relcoord[i] = new double[bodies*4];
	}
	//force provides fx,fy,fz for each body at step i and i+1
	double** force = new double*[2];
	for(i=0;i<2;i++){
		force[i] = new double[bodies*3];
	}
	//initialize matrix to zero
	//(note: gravityforces function does this each time it calls--no need to initialize)
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies*4;j++){
			relcoord[i][j]=0.0;
		}
	}
	//reinitialize relcoord to initial conditions
	rel_position(relcoord, output, 0, bodies);
	gravityforces(force, relcoord, mass, 0, bodies);

	//begin solution for-loop
	for(i=0;i<steps;i++){
		//Set next time value
		output[i+1][0] = output[i][0] + h;
		//calculate x,y,z		
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+1] = output[i][j*8+1] + h*output[i][j*8+5] + halfhsq*force[0][j*3];
			output[i+1][j*8+2] = output[i][j*8+2] + h*output[i][j*8+6] + halfhsq*force[0][j*3+1];
			output[i+1][j*8+3] = output[i][j*8+3] + h*output[i][j*8+7] + halfhsq*force[0][j*3+2];
			output[i+1][j*8+4] = sqrt(output[i+1][j*8+1]*output[i+1][j*8+1] + output[i+1][j*8+2]*output[i+1][j*8+2] + output[i+1][j*8+3]*output[i+1][j*8+3]);
		}
		//new relcoord
		rel_position(relcoord, output, i+1, bodies);
		gravityforces(force, relcoord, mass, 1, bodies);
		//calculate vx, vy, vz
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+5] = output[i][j*8+5] + halfh*(force[1][j*3] + force[0][j*3]);
			output[i+1][j*8+6] = output[i][j*8+6] + halfh*(force[1][j*3+1] + force[0][j*3+1]);
			output[i+1][j*8+7] = output[i][j*8+7] + halfh*(force[1][j*3+2] + force[0][j*3+2]);
			output[i+1][j*8+8] = sqrt(output[i+1][j*8+5]*output[i+1][j*8+5] + output[i+1][j*8+6]*output[i+1][j*8+6] + output[i+1][j*8+7]*output[i+1][j*8+7]);
		}
		//force_i for next loop is force_i+1 for this loop. No need for double computation.
		for(j=0;j<bodies;j++){
			force[0][j*3] = force[1][j*3];
			force[0][j*3+1] = force[1][j*3+1];
			force[0][j*3+2] = force[1][j*3+2];
		}
	}
	//delete dross--keep output
	for(i=0;i<bodies;i++){
		delete[] relcoord[i];
	}
	delete[] relcoord;
	for(i=0;i<2;i++){
		delete[] force[i];
	}
	delete[] force;

}
void verlet_relcor(double**& output, double*& mass, int bodies, int steps, double tmax, int sun_choice){
	int i, j, k, m, n;
	double h, halfh, halfhsq, rcubed, fourpisq;
	h = (tmax-0.0)/((double) (steps));
	halfh = 0.5*h;
	halfhsq = halfh*h;
	fourpisq = 4.0*M_PI*M_PI;
	//relcoord to provide relative coordinates, x, y, z, r, of body i with respect to body j
	//technically double calculating--consider revising
	double** relcoord = new double*[bodies];
	for(i=0;i<bodies;i++){
		relcoord[i] = new double[bodies*8];
	}
	double** relangmomentum = new double*[bodies];
	for(i=0;i<bodies;i++){
		relangmomentum[i] = new double[bodies*4];
	}
	//force provides fx,fy,fz for each body at step i and i+1
	double** force = new double*[2];
	for(i=0;i<2;i++){
		force[i] = new double[bodies*3];
	}
	//initialize matrices to zero
	//(note: gravityforces function does this each time it calls--no need to initialize)
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies*8;j++){
			relcoord[i][j]=0.0;
		}
	}
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies*4;j++){
			relangmomentum[i][j]=0.0;
		}
	}

	//initialize to initial conditions
	rel_coord(relcoord, output, 0, bodies);
	rel_angmomentum(relangmomentum, relcoord, bodies);
	gravityforces_relcorrection(force, relangmomentum, relcoord, mass, 0, bodies);

	//begin solution for-loop
	for(i=0;i<steps;i++){
		//Set next time value
		output[i+1][0] = output[i][0] + h;
		//calculate x,y,z		
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+1] = output[i][j*8+1] + h*output[i][j*8+5] + halfhsq*force[0][j*3];
			output[i+1][j*8+2] = output[i][j*8+2] + h*output[i][j*8+6] + halfhsq*force[0][j*3+1];
			output[i+1][j*8+3] = output[i][j*8+3] + h*output[i][j*8+7] + halfhsq*force[0][j*3+2];
			output[i+1][j*8+4] = sqrt(output[i+1][j*8+1]*output[i+1][j*8+1] + output[i+1][j*8+2]*output[i+1][j*8+2] + output[i+1][j*8+3]*output[i+1][j*8+3]);
		}

		//update forces
		rel_coord(relcoord, output, i+1, bodies);
		rel_angmomentum(relangmomentum, relcoord, bodies);
		gravityforces_relcorrection(force, relangmomentum, relcoord, mass, 1, bodies);
		
		//calculate vx, vy, vz
		for(j=sun_choice;j<bodies;j++){
			output[i+1][j*8+5] = output[i][j*8+5] + halfh*(force[1][j*3] + force[0][j*3]);
			output[i+1][j*8+6] = output[i][j*8+6] + halfh*(force[1][j*3+1] + force[0][j*3+1]);
			output[i+1][j*8+7] = output[i][j*8+7] + halfh*(force[1][j*3+2] + force[0][j*3+2]);
			output[i+1][j*8+8] = sqrt(output[i+1][j*8+5]*output[i+1][j*8+5] + output[i+1][j*8+6]*output[i+1][j*8+6] + output[i+1][j*8+7]*output[i+1][j*8+7]);
		}
		//force_i for next loop is force_i+1 for this loop. No need for double computation.
		for(j=0;j<bodies;j++){
			force[0][j*3] = force[1][j*3];
			force[0][j*3+1] = force[1][j*3+1];
			force[0][j*3+2] = force[1][j*3+2];
		}
	}
	//delete dross--keep output
	for(i=0;i<bodies;i++){
		delete[] relcoord[i];
	}
	delete[] relcoord;
	for(i=0;i<2;i++){
		delete[] force[i];
	}
	delete[] force;

}

/***************************
helionstates functions to 
retrieve aphelion & peri-
helion data
***************************/
void helionstates(double**& output, double**& aphelion, double**& perihelion, int bodies, int steps){
	int i, j;
	//Initialize with beginning values
	for(i=0;i<bodies;i++){
		aphelion[i][0] = output[0][i*8+4];
		perihelion[i][0] = aphelion[i][0];
		aphelion[i][1] = output[0][0];
		perihelion[i][1] = aphelion[i][1];
	}
	//sort through all radii and find the largest and smallest for each body
	for(i=1;i<steps+1;i++){
		for(j=0;j<bodies;j++){
			if(output[i][j*8+4]>aphelion[j][0]){
				aphelion[j][0] = output[i][j*8+4];
				aphelion[j][1] = output[i][0];
			}
			if(output[i][j*8+4]<perihelion[j][0]){
				perihelion[j][0] = output[i][j*8+4];
				perihelion[j][1] = output[i][0];
			}
		}
	}
}
void helionstates_dynamic(double**& aphelion, int& aphguess, double**& perihelion, int& periguess, double**& output, int steps, int body, int range, int dec_or_inc){
	int i, j, stepcount, aphcount, pericount, rposition, holder0, holder1;
	double maxtemp[2], mintemp[2], holder[2];
	rposition = body*8+4;

	for(i=1;i<5;i++){
		aphelion[0][i] = output[0][body*8+i];	
		perihelion[0][i] = aphelion[0][i];
	}

	aphelion[0][0] = output[0][0];
	perihelion[0][0] = aphelion[0][0];
	holder[0] = aphelion[0][4];
	holder[1] = 0;

	stepcount = 0;
	aphcount=1;
	pericount=1;	
	
	while(stepcount<steps+1){
		if(aphcount >= aphguess){
			matrix_resize(aphelion, aphguess, 5, aphcount+1, 5);
			aphelion[aphcount][0] = -1;
			for(i=0;i<5;i++){
				aphelion[aphcount][i] = 0;
			}
			aphguess = aphcount+1;
		}if(pericount >= periguess){
			matrix_resize(perihelion, periguess, 5, pericount+1, 5);
			perihelion[pericount][0] = -1;
			for(i=1;i<5;i++){
				perihelion[pericount][i] = 0;
			}
			periguess = pericount+1;
		}

		if(dec_or_inc==0){
			mintemp[0] = holder[0];
			mintemp[1] = holder[1];
			for(i=stepcount;i<stepcount+range;i++){
				if(output[i][rposition]<=mintemp[0]){
					mintemp[0] = output[i][rposition];
					mintemp[1] = i;
				}
			}
			if(mintemp[0]<holder[0]){
				holder[0] = mintemp[0];
				holder[1] = mintemp[1];
			}
			else{
				dec_or_inc = 1;
				holder1=holder[1];
				perihelion[pericount][0] = output[holder1][0];
				for(i=1;i<5;i++){
					perihelion[pericount][i] = output[holder1][body*8+i];
				}	
				pericount++;
			}
		}
		else if(dec_or_inc==1){
			maxtemp[0] = holder[0];
			mintemp[1] = holder[1];
			for(i=stepcount;i<stepcount+range;i++){
				if(output[i][rposition]>=maxtemp[0]){
					maxtemp[0] = output[i][rposition];
					maxtemp[1] = i;
				}
			}
			if(maxtemp[0]>holder[0]){
				holder[0] = maxtemp[0];
				holder[1] = maxtemp[1];
			}
			else{
				dec_or_inc = 0;
				holder1=holder[1];
				aphelion[aphcount][0] = output[holder1][0];
				for(i=1;i<5;i++){			
					aphelion[aphcount][i] = output[holder1][body*8+i];
				}
				aphcount++;
			}
		}
		stepcount += range;
		if(stepcount + range > steps+1){
			range = (steps+1) - stepcount;
		}
	}
}

#endif
