/****************************************************
This code file solves the evolution of a system 
of massive bodies given their respective mass, 
positions, and velocities. In particular, we 
model the solar system given data from Nasa's Jet
Propulsion Laboratory. We proceed by two methods, 
Verlet and Runge Kutta to the fourth order. 


Verlet takes advantage of Taylor expansions, 
resulting in the following equations for the
positions and velocities:

x_i+1 = x_i + h*v_i + (h^2)/2*Force_i/mass
v_i+1 = v_i + (h/2)*(Force_i+1 + Force_i)/mass


RK4 uses Simpson's rule, that is
integral_(t)^(t+h) (f(x,t)) 
 = (h/6)*(f(i) + 4*f(i+1/2) + f(i+1))
 = (h/6)*(f(i) + 2*f(i+1/2) + 2*f(i+1/2) + f(i+1)
 = (h/6)*(k1 + 2*k2 + 2*k3 + k4),
where we define f(i) = f(x_i, t_i). In our case, 
we must perform this operation for both x and v, 
where we follow
 
	k1v = force(i)/mass
	k1x = v(i)
	x_i+1/2 = x(i) + (h/2)*k1x

	k2v = force(i+1/2)/mass
	k2x = v(i) + (h/2)*k2v
	x_i+1/2 = x(i) + (h/2)*k2x

	k3v = force(i+1/2)/mass
	k3x = v(i) + (h/2)*k3v
	x_i+1 = x(i) + h*k3x

	k4v = force(i+1)/mass
	k4x = v(i) + h*k4v

	v_i+1 = v(i) + (h/6)*(k1v + 2*k2v + 2*k3v + k4v)
	x_i+1 = x(i) + (h/6)*(k1x + 2*k2x + 2*k3x + k4x)


There is also the option to allow the system to
evolve using the first relativistic correction to the
force. The methods follow exactly as above once this
modification has been made to the force.

****************************************************/

#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include "time.h"
#include "project3_library.h"

using namespace std;

ofstream ofile, ofile2;

int main(int argc, char* argv[]){
	double** output, **aphelion, **perihelion;
	char* outfilename;

	int i, j, k, steps, planet_choice, sun_choice, method, outputtype, bodies;
	int stepsperyear, timesperyear;
	int dec_or_inc, aphguess, periguess;
	double tmax, time, timesort;

	clock_t start, finish;
	string periname, aphname;
	stringstream number;

	if(argc<8){
		cout << "Bad usage. Enter also 'outfilename planet_choice immobilize_sun_choice method_choice output_type steps tmax' on same line." << endl;
		cout << "Current options for planet_choice include" << endl;
		cout << "\t1: coplanar earth-sun system" << endl << "\t2: coplanar earth-jupiter-sun system" << endl;
		cout << "\t3: full solar system, with planets and luna" << endl;
		cout << "\t4: coplanar sun-mercury system" << endl << "\t5: rogue star full solar system, for fun" << endl;
		cout << "Current options for immobilze_sun_choice include" << endl;
		cout << "\t0: allow sun to evolve in coordinate space" << endl << "\t1: sun is initialized only and is not affected by other bodies" << endl;
		cout << "Current options for method_choice include" << endl;
		cout << "\t1: RK4" << endl << "\t2: Verlet" << endl << "\t3: non-class Verlet" << endl << "\t4: non-class relativistic Verlet" << endl;
		cout << "Current options for output_type include" << endl;
		cout << "\t1: pure output, steps by 8*bodies + 1" << endl << "\t2: as one, with print maximum aphelion and minimum perihelion to screen" << endl;
		cout << "\t3: 'timesperyear' ouput per earth year, enter after 'tmax'" << endl << "\t4: as three, with helionstates data files created" << endl;
		exit(1);
	}
	else{
		outfilename = argv[1];
		planet_choice = atoi(argv[2]);
		sun_choice = atoi(argv[3]);
		method = atoi(argv[4]);
		outputtype = atoi(argv[5]);
		steps = atoi(argv[6]);
		tmax = atof(argv[7]);
		if(outputtype==3||outputtype==4){
			if(argc<9){
				cout << "Bad usage. Enter also 'timesperyear' after 'tmax.'" << endl;
				exit(1);
			}
			tmax = atoi(argv[7]);
			timesperyear = atoi(argv[8]);
			if(timesperyear==0){
				cout << "Bad usage. Need a non-zero 'timesperyear.'" << endl;
				exit(1);
			}
			if(timesperyear*tmax>steps){
				cout << "Bad usage. 'timesperyear'*'tmax' is bigger than total 'steps.' Insufficient data points for output." << endl;
				exit(1);
			}
			stepsperyear = steps/(tmax*timesperyear);
		}	
	}
		
/********************************************
Here are the major massive bodies present 
within our solar system. Comment out any
which are not desired for the calculation.

Current options are
1: coplanar earth-sun system
2: coplanar earth-jupiter-sun system
3: full solar system, with planets and luna
4: coplanar sun-mercury system
5: rogue star full solar system, for fun
********************************************/

		massivesystem solarsystem;

	if(planet_choice==1){
		massivebody sun(1,0,0,0,0,0,0);
		massivebody earth(3.003489e-6,1,0,0,0,2.0*M_PI,0);
		solarsystem.add(sun);
		solarsystem.add(earth);
	}
	if(planet_choice==2){
		if(sun_choice==0){
			massivebody sun(1,-3.003489e-6,9.5479194e-4*5.2,0,-9.5479194e-4*0.439*2.0*M_PI,-3.003489e-6*2.0*M_PI,0);
			massivebody earth(3.003489e-6,1,0,0,0,2.0*M_PI,0);
			massivebody jupiter(9.5479194e-4,0,-5.2,0,0.439*2.0*M_PI,0,0);
			solarsystem.add(sun);
			solarsystem.add(earth);
			solarsystem.add(jupiter);
		}		
		if(sun_choice==1){
			massivebody sun(1,0,0,0,0,0,0);
			massivebody earth(3.003489e-6,1,0,0,0,2.0*M_PI,0);
			massivebody jupiter(9.5479194e-4,0,-5.2,0,0.439*2.0*M_PI,0,0);
			solarsystem.add(sun);
			solarsystem.add(earth);
			solarsystem.add(jupiter);
		}
	}
	if(planet_choice==3){
		massivebody sun(1,3.771551748320805E-03, 1.938413234187417E-03, -1.625928558000791E-04,365*3.428095941073785E-08, 365*6.978886434512168E-06, 365*-9.372671992938156E-09);
		massivebody mercury(1.6601e-7,3.566221110752382E-01,-1.449153604767920E-01,-4.453344798939488E-02, 365*5.324168987021533E-03,365*2.725352157519689E-02, 365*1.737877336062238E-03);
		massivebody venus(2.4478383e-6,4.456189829066335E-01,-5.759198980424926E-01,-3.358283335789598E-02,365*1.593371764525070E-02, 365*1.222110108236829E-02,365*-7.520174121955455E-04);
		massivebody earth(3.003489e-6,-9.906650404586314E-01,4.353612431574581E-02,-1.569714899841466E-04,365*-9.933032846586867E-04,365*-1.724423592380582E-02,365*2.880574493748607E-07);
		massivebody moon(1.230004e-8, -9.917828800119005E-01, 4.587463373520608E-02, -3.563694480647205E-04, 365*-1.530892843498182E-03, 365*-1.746686609878291E-02, 365*2.803200805250705E-05);
		massivebody mars(3.227151e-7, -1.394909885799664E+00, -7.759974369033253E-01, 1.786251125355763E-02, 365*7.324673346484570E-03, 365*-1.102624283521118E-02, 365*-4.109846883854566E-04);
		massivebody jupiter(9.5479194e-4,-5.321136962878863E+00, 1.055810040982731E+00, 1.146123452783326E-01, 365*-1.556597232697437E-03, 365*-7.046207863842619E-03, 365*6.409351264039102E-05);
		massivebody saturn(2.858860e-4,-3.332098484988519E+00, -9.442038663142483E+00, 2.967846224795282E-01, 365*4.954645896896637E-03, 365*-1.873255191158998E-03, 365*-1.647257228743735E-04);
		massivebody uranus(4.366244e-5, 1.876841090026478E+01, 6.826065082612249E+00, -2.177966843356428E-01, 365*-1.372937790621792E-03, 365*3.512867479374193E-03, 365*3.086802162915850E-05);
		massivebody neptune(5.151389e-5, 2.803900494548452E+01, -1.054870089186826E+01, -4.289565554838171E-01, 365*1.084650993604757E-03, 365*2.957157649376530E-03, 365*-8.562727609311126E-05);
		massivebody pluto(6.5812e-9, 8.773368933896196E+00, -3.186331328356860E+01, 8.718065633574812E-01, 365*3.100891963853092E-03, 365*1.939401372093854E-04, 365*-9.194995916567601E-04);

		solarsystem.add(sun);
		solarsystem.add(mercury);
		solarsystem.add(venus);
		solarsystem.add(earth);
		solarsystem.add(moon);
		solarsystem.add(mars);
		solarsystem.add(jupiter);
		solarsystem.add(saturn);
		solarsystem.add(uranus);
		solarsystem.add(neptune);
		solarsystem.add(pluto);
	}
	if(planet_choice==4){
		massivebody sun(1,0,0,0,0,0,0);
		massivebody mercury(1.6601e-7,0.3075, 0, 0, 0, 12.44, 0);
	
		solarsystem.add(sun);
		solarsystem.add(mercury);
	}
	if(planet_choice==5){

		massivebody sun(1,3.771551748320805E-03, 1.938413234187417E-03, -1.625928558000791E-04,365*3.428095941073785E-08, 365*6.978886434512168E-06, 365*-9.372671992938156E-09);
		massivebody mercury(1.6601e-7,3.566221110752382E-01,-1.449153604767920E-01,-4.453344798939488E-02, 365*5.324168987021533E-03,365*2.725352157519689E-02, 365*1.737877336062238E-03);
		massivebody venus(2.4478383e-6,4.456189829066335E-01,-5.759198980424926E-01,-3.358283335789598E-02,365*1.593371764525070E-02, 365*1.222110108236829E-02,365*-7.520174121955455E-04);
		massivebody earth(3.003489e-6,-9.906650404586314E-01,4.353612431574581E-02,-1.569714899841466E-04,365*-9.933032846586867E-04,365*-1.724423592380582E-02,365*2.880574493748607E-07);
		massivebody moon(1.230004e-8, -9.917828800119005E-01, 4.587463373520608E-02, -3.563694480647205E-04, 365*-1.530892843498182E-03, 365*-1.746686609878291E-02, 365*2.803200805250705E-05);
		massivebody mars(3.227151e-7, -1.394909885799664E+00, -7.759974369033253E-01, 1.786251125355763E-02, 365*7.324673346484570E-03, 365*-1.102624283521118E-02, 365*-4.109846883854566E-04);
		massivebody jupiter(9.5479194e-4,-5.321136962878863E+00, 1.055810040982731E+00, 1.146123452783326E-01, 365*-1.556597232697437E-03, 365*-7.046207863842619E-03, 365*6.409351264039102E-05);
		massivebody saturn(2.858860e-4,-3.332098484988519E+00, -9.442038663142483E+00, 2.967846224795282E-01, 365*4.954645896896637E-03, 365*-1.873255191158998E-03, 365*-1.647257228743735E-04);
		massivebody uranus(4.366244e-5, 1.876841090026478E+01, 6.826065082612249E+00, -2.177966843356428E-01, 365*-1.372937790621792E-03, 365*3.512867479374193E-03, 365*3.086802162915850E-05);
		massivebody neptune(5.151389e-5, 2.803900494548452E+01, -1.054870089186826E+01, -4.289565554838171E-01, 365*1.084650993604757E-03, 365*2.957157649376530E-03, 365*-8.562727609311126E-05);
		massivebody pluto(6.5812e-9, 8.773368933896196E+00, -3.186331328356860E+01, 8.718065633574812E-01, 365*3.100891963853092E-03, 365*1.939401372093854E-04, 365*-9.194995916567601E-04);
		massivebody roguestar(1,-500,0,0,50,0,0);

		solarsystem.add(sun);
		solarsystem.add(mercury);
		solarsystem.add(venus);
		solarsystem.add(earth);
		solarsystem.add(moon);
		solarsystem.add(mars);
		solarsystem.add(jupiter);
		solarsystem.add(saturn);
		solarsystem.add(uranus);
		solarsystem.add(neptune);
		solarsystem.add(pluto);
		solarsystem.add(roguestar);
	}

	bodies = solarsystem.massivebody_count;
	matrix_alloc(output, steps+1, bodies*8+1);

/**********************************************
method choice
	1: RK4
	2: Verlet
	3: Non-massivesystem Verlet
	4: Verlet with relativistic correction
**********************************************/
	if(method==1){
		start = clock();
		solarsystem.RK4(output,steps,tmax,sun_choice);
		finish = clock();
	}
	else if(method==2){
		start = clock();
		solarsystem.verlet(output,steps,tmax,sun_choice);
		finish = clock();
	}
	else if(method==3){
		double* mass = new double[bodies];
		start = clock();
		solarsystem.initialize(output, mass, steps);
		verlet(output,mass,bodies,steps,tmax,sun_choice);
		finish = clock();
		delete[] mass;
	}
	else if(method==4){
		double* mass = new double[bodies];
		start = clock();
		solarsystem.initialize(output, mass, steps);
		verlet_relcor(output,mass,bodies,steps,tmax,sun_choice);
		finish = clock();
		delete[] mass;
	}
	time = (finish - start)/((double) CLOCKS_PER_SEC);

/*****************************************************************
outputfile choice
	1: Pure output file. No modifications
	2: As one, but with the largest and smallest, respectively, 
		aphelion and perihelion for each body printed to screen
	3: Output 'timesperyear' data points per year, for long 
		time spans
	4: As two, but with data files for each body containing all 
		aphelion and perilion positions
*****************************************************************/

	ofile.open(outfilename);
	ofile.precision(8);
	if(outputtype==1){
		start = clock();
		for(i=0;i<steps+1;i++){
			for(j=0;j<bodies*8+1;j++){
				ofile << setw(15) << output[i][j];
			}
			ofile << endl;
		}
		finish = clock();
	}
	else if(outputtype==2){
		start = clock();
		matrix_alloc(aphelion,bodies,2);
		matrix_alloc(perihelion,bodies,2);
		for(i=0;i<steps+1;i++){
			for(j=0;j<bodies*8+1;j++){
				ofile << setw(15) << output[i][j];
			}
			ofile << endl;
		}
		helionstates(output, aphelion, perihelion, bodies, steps);
		for(i=0;i<bodies;i++){
				cout << setw(15) << aphelion[i][0] << setw(15) << aphelion[i][1];
				cout << setw(30) << perihelion[i][0] << setw(15) << perihelion[i][1];
				cout << endl;
		}
		matrix_delete(aphelion,bodies);
		matrix_delete(perihelion,bodies);
		finish = clock();
	}
	else if(outputtype==3){
		start = clock();
		for(i=0;i<steps+1;i+=stepsperyear){
			for(j=0;j<bodies*8+1;j++){
				ofile << setw(15) << output[i][j];
			}
			ofile << endl;
		}
		finish = clock();
	}
	else if(outputtype==4){
		start = clock();
		for(i=0;i<steps+1;i+=stepsperyear){
			for(j=0;j<bodies*8+1;j++){
				ofile << setw(15) << output[i][j];
			}
			ofile << endl;
		}
		for(k=0;k<bodies;k++){
			periname = "perihelion_body_";
			aphname = "aphelion_body_";
			aphguess = 1;
			periguess = aphguess;
			number << k;
			periname += number.str() + ".dat";
			aphname += number.str() + ".dat";
			number.str("");
			matrix_alloc(aphelion,aphguess,5);
			matrix_alloc(perihelion,periguess,5);
			dec_or_inc=1;
			helionstates_dynamic(aphelion, aphguess, perihelion, periguess, output,steps,k,20,dec_or_inc);
			
			ofile2.open(periname);
			ofile2.precision(8);
			for(i=0;i<periguess;i++){
				for(j=0;j<5;j++){
					ofile2 << setw(15) << perihelion[i][j];
				}
				ofile2 << endl;
			}
			ofile2.close();
			ofile2.open(aphname);
			ofile2.precision(8);
			for(i=0;i<aphguess;i++){
				for(j=0;j<5;j++){
					ofile2 << setw(15) << aphelion[i][j];
				}
				ofile2 << endl;
			}
			ofile2.close();

			matrix_delete(aphelion,aphguess);
			matrix_delete(perihelion,periguess);
		}
		finish = clock();
	}
	timesort = (finish - start)/((double) CLOCKS_PER_SEC);

	//print time for calculation to screen
	cout << "\tcomputation time " << time << " seconds" << endl;
	cout << "\toutput time " << timesort << " seconds" << endl;
	ofile.close();

	matrix_delete(output,steps+1);
	return 0;
}
