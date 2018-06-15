#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <limits> 
#include <iomanip> 
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <string.h>
//#include <cstdlib>
using namespace std;
string preference;
const int n_c = 6, number_of_particles = pow(n_c,2), u=sqrt(number_of_particles); 
double  a ,h, L_x=7, L_y=sqrt(3)*L_x/2;//L_x, L_y;
const int step = 1000;//skip 1000 coordinates
const int NTIME=1000000 /* number of steps */, control_temp=0 /* T control enable/disable 1/0 */;
int o = 0; //T control steps
double length[number_of_particles][number_of_particles];
const double epsilon=0.0031, lattice_constant=40, dt=0.001, primary_velocity = 1;
double temperature_average, temperature_desired, temperature_prev, temperature, temperature_sum, kin_energy[NTIME], pot_energy[NTIME];
int stop = 0;
double R, Z, dx, dy, sigma=0.01;
const int control_point_1=10/dt, control_point_2=50/dt, control_point_3=100/dt; //steps at which temperature control starts
double random(double min, double max)
{
    return (double)(rand())/RAND_MAX*(max - min) + min;
}
//STRUCTURE
struct primary_pairs
{
	public:
		double x;
		double y;
};
 	primary_pairs coord[number_of_particles][NTIME], acceleration[number_of_particles][NTIME], velocity[number_of_particles][NTIME];
	primary_pairs delta[number_of_particles][number_of_particles], velocity_average, sum_of_velocities;
int main() 
{

ofstream pulse, write, kinwrite, potwrite, fullwrite, tempwrite, avtempwrite, partwrite, rwrite, zwrite;
ofstream triangle5, triangle7, square5, square7;
ofstream write_tr_5, write_tr_7, write_sq_5, write_sq_7;
write_tr_5.open("write_tr_5");
write_tr_7.open("write_tr_7");
write_sq_5.open("write_sq_5");
write_sq_7.open("write_sq_7");
write.open("data.txt");
kinwrite.open("kinenergy");
potwrite.open("potenergy");
fullwrite.open("fullenergy");
tempwrite.open("temperature");
avtempwrite.open("avtemperature");
partwrite.open("particle");
rwrite.open("R");
zwrite.open("Z");
pulse.open("pulse");
triangle5.open("triangle5");
triangle7.open("triangle7");
square5.open("square5");
square7.open("square7");
		preference="square5";
		L_x=5;
		L_y=L_x*(sqrt(3)/2);
	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			length[i][j]=0;
		}
	}	
	//RANDOM NUMBER GENERATOR
    srand((unsigned int)time(0));

	for (int i=0; i<number_of_particles; i++)
	{
		for (int N=0; N<NTIME; N++)
		{
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			acceleration[i][N].x = 0;
			acceleration[i][N].y = 0;
			velocity[i][N].x = 0;
			velocity[i][N].y = 0;
			kin_energy[N] = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}

	
//COORDINATES 0
//for (int i=0; i<4; i++){
	
///////////////
	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			length[i][j]=0;
		}
	}	

	for (int i=0; i<number_of_particles; i++)
	{
		for (int N=0; N<NTIME; N++)
		{
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			acceleration[i][N].x = 0;
			acceleration[i][N].y = 0;
			velocity[i][N].x = 0;
			velocity[i][N].y = 0;
			kin_energy[N] = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}
////////////////
	double sigma_min=0, min_energy=10000;
	pot_energy[0]=0;
	sigma=0;

	cout << preference << ", L_x = " << L_x << endl;

	if(preference=="square5" or preference=="square7"){
		a=sqrt(L_x*L_y)/(n_c);
		for (int i=0; i<u;  i++)
		{
			for (int j=0; j<u; j++){
				coord[i*u+j][0].x= j*a;
				coord[i*u+j][0].y=i*a;
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
				if(preference=="square5"){
				write_sq_5 <<  coord[i*u+j][0].x << " " << coord[i*u+j][0].y << endl;
				}
				if(preference=="square7"){
				write_sq_7 <<  coord[i*u+j][0].x << " " << coord[i*u+j][0].y << endl;
				}
			}
		}

	}
		//cout << "" << endl;
		//cout << "" << endl;
		write << "" << endl;
		write << "" << endl;
		for (int i=0; i<number_of_particles; i++) //CHOOSE A PARTICLE TO FIND A DISTANCE
		{
				for ( int j=0; j<number_of_particles; j++)
				{
					
					if (i==j)
					{
						continue;
					}
					delta[i][j].x=(coord[i][0].x-coord[j][0].x);
					delta[i][j].y=(coord[i][0].y-coord[j][0].y);
					
					if (abs(delta[i][j].x)>L_x/2 )
						{
							delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-L_x);
						}
					if (abs(delta[i][j].y)>L_y/2 )
						{
							delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-L_y);
						}
					length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
				}
		}
		while (sigma < 5){
			for (int i=0; i<number_of_particles; i++) {
				for ( int j=i+1; j<number_of_particles; j++){
					pot_energy[0]+=4*epsilon*( pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
				}
			}
			//cout <<  "potential energy, sigma = " << sigma << " : " <<  pot_energy[0] << setprecision(15)  << endl;
			square5 << sigma << " " << pot_energy[0] << setprecision(15)  << endl;
			if (pot_energy[0]<min_energy){
				min_energy=pot_energy[0];
				sigma_min=sigma;
				}
			
			pot_energy[0]=0;
			sigma+=0.001;
		}
		cout << "MIN_ENERGY = " << min_energy << endl; 
		cout << "SIGMA = " << sigma_min << endl; 
		cout << "DENSITY = " << pow(n_c,2)/(L_x*L_y);
		cout << endl;
		pot_energy[0]=min_energy;
		sigma=sigma_min;
		//sigma=0.74;
	sigma_min=0;
	min_energy=10000;
	pot_energy[0]=0;
	sigma=0;
		cout << endl;
		///////////////////////////////////////
		///////////////////////////////////////		
		///////////////////////////////////////		
		
		
		preference="square7";
		L_x=7;
		L_y=L_x*(sqrt(3)/2);



	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			length[i][j]=0;
		}
	}	
	//RANDOM NUMBER GENERATOR
    srand((unsigned int)time(0));

	for (int i=0; i<number_of_particles; i++)
	{
		for (int N=0; N<NTIME; N++)
		{
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			acceleration[i][N].x = 0;
			acceleration[i][N].y = 0;
			velocity[i][N].x = 0;
			velocity[i][N].y = 0;
			kin_energy[N] = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}
  	

//COORDINATES 0
//for (int i=0; i<4; i++){
	
///////////////
	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			length[i][j]=0;
		}
	}	

	for (int i=0; i<number_of_particles; i++)
	{
		for (int N=0; N<NTIME; N++)
		{
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			acceleration[i][N].x = 0;
			acceleration[i][N].y = 0;
			velocity[i][N].x = 0;
			velocity[i][N].y = 0;
			kin_energy[N] = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}
////////////////
	sigma_min=0;
	min_energy=10000;
	pot_energy[0]=0;
	sigma=0;
	cout << preference << ", L_x = " << L_x << endl;
	if(preference=="square5" or preference=="square7"){
		a=sqrt(L_x*L_y)/(n_c); 
		for (int i=0; i<u;  i++)
		{
			for (int j=0; j<u; j++)
			{
				coord[i*u+j][0].x= j*a;
				coord[i*u+j][0].y=i*a;
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
				if(preference=="square7"){
				write_sq_7 <<  coord[i*u+j][0].x << " " << coord[i*u+j][0].y << endl;
			}
			}
		}

	}
		//cout << "" << endl;
		//cout << "" << endl;
		write << "" << endl;
		write << "" << endl;
		for (int i=0; i<number_of_particles; i++) //CHOOSE A PARTICLE TO FIND A DISTANCE
		{
				for ( int j=0; j<number_of_particles; j++)
				{
					
					if (i==j)
					{
						continue;
					}
					delta[i][j].x=(coord[i][0].x-coord[j][0].x);
					delta[i][j].y=(coord[i][0].y-coord[j][0].y);
					
					if (abs(delta[i][j].x)>L_x/2 )
						{
							delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-L_x);
						}
					if (abs(delta[i][j].y)>L_y/2 )
						{
							delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-L_y);
						}
					length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
				}
		}
		while (sigma < 5){
			for (int i=0; i<number_of_particles; i++) {
				for ( int j=i+1; j<number_of_particles; j++){
					pot_energy[0]+=4*epsilon*( pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
				}
			}
			//cout <<  "potential energy, sigma = " << sigma << " : " <<  pot_energy[0] << setprecision(15)  << endl;
			square7 << sigma << " " << pot_energy[0] << setprecision(15)  << endl;
			if (pot_energy[0]<min_energy){
				min_energy=pot_energy[0];
				sigma_min=sigma;
				}
			
			pot_energy[0]=0;
			sigma+=0.001;
		}
		cout << "MIN_ENERGY = " << min_energy << endl; 
		cout << "SIGMA = " << sigma_min << endl;
		cout << "DENSITY = " << pow(n_c,2)/(L_x*L_y);
		cout << endl;
		pot_energy[0]=min_energy;
		sigma=sigma_min;
		//sigma=0.74;
		cout << endl;
		
		
		
		
		
		///////////////////////////////////////
		///////////////////////////////////////		
		///////////////////////////////////////		
		
		
		preference="triangle5";
		L_x=5;
		L_y=L_x*(sqrt(3)/2);
	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			length[i][j]=0;
		}
	}	
	//RANDOM NUMBER GENERATOR
    srand((unsigned int)time(0));

	for (int i=0; i<number_of_particles; i++)
	{
		for (int N=0; N<NTIME; N++)
		{
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			acceleration[i][N].x = 0;
			acceleration[i][N].y = 0;
			velocity[i][N].x = 0;
			velocity[i][N].y = 0;
			kin_energy[N] = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}
//COORDINATES 0
//for (int i=0; i<4; i++){
	
///////////////
	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			length[i][j]=0;
		}
	}	

	for (int i=0; i<number_of_particles; i++)
	{
		for (int N=0; N<NTIME; N++)
		{
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			acceleration[i][N].x = 0;
			acceleration[i][N].y = 0;
			velocity[i][N].x = 0;
			velocity[i][N].y = 0;
			kin_energy[N] = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}
////////////////
	sigma_min=0;
	min_energy=10000;
	pot_energy[0]=0;
	sigma=0;
	cout << preference << ", L_x = " << L_x << endl;
	if (preference=="triangle5" or preference=="triangle7"){
		a=L_x/(n_c);
		h = a*(sqrt(3)/2);
		for (int i=0; i<u;  i++)
		{
			for (int j=0; j<u; j++)
			{
				coord[i*u+j][0].x= i%2==0  ? j*a : j*a+a/2;
				coord[i*u+j][0].y=i*h;
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
				if(preference=="triangle5"){
				write_tr_5 <<  coord[i*u+j][0].x << " " << coord[i*u+j][0].y << endl;
			}
			}
		}

	}
		//cout << "" << endl;
		//cout << "" << endl;
		write << "" << endl;
		write << "" << endl;
		for (int i=0; i<number_of_particles; i++) //CHOOSE A PARTICLE TO FIND A DISTANCE
		{
				for ( int j=0; j<number_of_particles; j++)
				{
					
					if (i==j)
					{
						continue;
					}
					delta[i][j].x=(coord[i][0].x-coord[j][0].x);
					delta[i][j].y=(coord[i][0].y-coord[j][0].y);
					
					if (abs(delta[i][j].x)>L_x/2 )
						{
							delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-L_x);
						}
					if (abs(delta[i][j].y)>L_y/2 )
						{
							delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-L_y);
						}
					length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
				}
		}
		while (sigma < 5){
			for (int i=0; i<number_of_particles; i++) {
				for ( int j=i+1; j<number_of_particles; j++){
					pot_energy[0]+=4*epsilon*( pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
				}
			}
			//cout <<  "potential energy, sigma = " << sigma << " : " <<  pot_energy[0] << setprecision(15)  << endl;
			if (preference=="triangle5"){
				triangle5 << sigma << " " << pot_energy[0] << setprecision(15)  << endl;
				}
			if (pot_energy[0]<min_energy){
				min_energy=pot_energy[0];
				sigma_min=sigma;
				}
			
			pot_energy[0]=0;
			sigma+=0.001;
		}
		cout << "MIN_ENERGY = " << min_energy << endl; 
		cout << "SIGMA = " << sigma_min << endl;
		cout << "DENSITY = " << pow(n_c,2)/(L_x*L_y);
		cout << endl;
		pot_energy[0]=min_energy;
		sigma=sigma_min;
	sigma_min=0;
	min_energy=10000;
	pot_energy[0]=0;
	sigma=0;
		//sigma=0.74;
		cout << endl;

		///////////////////////////////////////
		///////////////////////////////////////		
		///////////////////////////////////////		
		
		
		preference="triangle7";
		L_x=7;
		L_y=L_x*(sqrt(3)/2);
		
	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			length[i][j]=0;
		}
	}	
	//RANDOM NUMBER GENERATOR
    srand((unsigned int)time(0));

	for (int i=0; i<number_of_particles; i++)
	{
		for (int N=0; N<NTIME; N++)
		{
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			acceleration[i][N].x = 0;
			acceleration[i][N].y = 0;
			velocity[i][N].x = 0;
			velocity[i][N].y = 0;
			kin_energy[N] = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}
 //COORDINATES 0
//for (int i=0; i<4; i++){
	
///////////////
	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			length[i][j]=0;
		}
	}	

	for (int i=0; i<number_of_particles; i++)
	{
		for (int N=0; N<NTIME; N++)
		{
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			acceleration[i][N].x = 0;
			acceleration[i][N].y = 0;
			velocity[i][N].x = 0;
			velocity[i][N].y = 0;
			kin_energy[N] = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}
////////////////
	sigma_min=0;
	min_energy=10000;
	pot_energy[0]=0;
	sigma=0;
	cout << preference << ", L_x = " << L_x << endl;
	if (preference=="triangle5" or preference=="triangle7"){
		a=L_x/(n_c);
		h = a*(sqrt(3)/2); 
		for (int i=0; i<u;  i++)
		{
			for (int j=0; j<u; j++)
			{
				coord[i*u+j][0].x= i%2==0  ? j*a : j*a+a/2;
				coord[i*u+j][0].y=i*h;
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
			write_tr_7 <<  coord[i*u+j][0].x << " " << coord[i*u+j][0].y << endl;
			}
		}

	}
		//cout << "" << endl;
		//cout << "" << endl;
		write << "" << endl;
		write << "" << endl;
		for (int i=0; i<number_of_particles; i++) //CHOOSE A PARTICLE TO FIND A DISTANCE
		{
				for ( int j=0; j<number_of_particles; j++)
				{
					
					if (i==j)
					{
						continue;
					}
					delta[i][j].x=(coord[i][0].x-coord[j][0].x);
					delta[i][j].y=(coord[i][0].y-coord[j][0].y);
					
					if (abs(delta[i][j].x)>L_x/2 )
						{
							delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-L_x);
						}
					if (abs(delta[i][j].y)>L_y/2 )
						{
							delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-L_y);
						}
					length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
				}
		}
		while (sigma < 5){
			for (int i=0; i<number_of_particles; i++) {
				for ( int j=i+1; j<number_of_particles; j++){
					pot_energy[0]+=4*epsilon*( pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
				}
			}
			//cout <<  "potential energy, sigma = " << sigma << " : " <<  pot_energy[0] << setprecision(15)  << endl;
			if (preference=="triangle7"){
				triangle7 << sigma << " " << pot_energy[0] << setprecision(15)  << endl;
				}
			if (pot_energy[0]<min_energy){
				min_energy=pot_energy[0];
				sigma_min=sigma;
				}
			
			pot_energy[0]=0;
			sigma+=0.001;
		}
		cout << "MIN_ENERGY = " << min_energy << endl; 
		cout << "SIGMA = " << sigma_min << endl;
		cout << "DENSITY = " << pow(n_c,2)/(L_x*L_y);
		cout << endl;
		pot_energy[0]=min_energy;
		sigma=sigma_min;
		//sigma=0.74;
	sigma_min=0;
	min_energy=10000;
	pot_energy[0]=0;
	sigma=0;
		cout << endl;
		
	
	//}


		//cout << "" << endl;
		//cout << "" << endl;

write.close();
kinwrite.close();
potwrite.close();
fullwrite.close();
tempwrite.close();
avtempwrite.close();
partwrite.close();
rwrite.close();
	return 0;
}
