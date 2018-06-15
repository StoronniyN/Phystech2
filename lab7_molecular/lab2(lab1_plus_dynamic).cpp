#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <limits> 
#include <iomanip> 
#include <math.h>
#include <fstream>
#include <stdio.h>
//#include <cstdlib>
using namespace std;

const int n_c = 6, number_of_particles = pow(n_c,2), u=sqrt(number_of_particles); 
double  L_x=7, L_y=L_x*(sqrt(3)/2), a=L_x/(n_c), h = a*(sqrt(3)/2); 
const int step = 1000; //skip 1000 coordinates
const int NTIME=1000000 /* number of steps */, control_temp=0 /* T control enable/disable 1/0 */;
int o = 0; //T control steps
double length[number_of_particles][number_of_particles];
const double epsilon=0.0031, lattice_constant=40, dt=0.0001, primary_velocity = 0.001, sigma=0.848351;
double temperature_average, temperature_desired, temperature_prev, temperature, temperature_sum, kin_energy[NTIME], pot_energy[NTIME];
int stop = 0;
double R, Z, dx, dy;
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
int main(){
ofstream pulse, write, kinwrite, potwrite, fullwrite, tempwrite, avtempwrite, partwrite, rwrite, zwrite, sigma_p_energy;
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
sigma_p_energy.open("p_energy");
	for (int i=0; i<number_of_particles; i++){
		for (int j=0; j<number_of_particles; j++){
			length[i][j]=0;
		}
	}	
	//RANDOM NUMBER GENERATOR
    srand((unsigned int)time(0));

	for (int i=0; i<number_of_particles; i++){
		for (int N=0; N<NTIME; N++){
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
    for (int i=0; i<number_of_particles; i++)
    {
		velocity[i][0].x =random(-primary_velocity, primary_velocity);
		//////cout << "v_x_" << i << "= " << velocity[i].x << endl;
		velocity[i][0].y = random(-primary_velocity, primary_velocity);
		//////cout << "v_y_" <<i << "= " << velocity[i].y << endl;
	}	
// SUM[v_i] - number_of_particles*<V> = 0, <V>=SUM[v_i]

	velocity_average.x = 0;
	velocity_average.y = 0;
	sum_of_velocities.x = 0;
	sum_of_velocities.y = 0;
	for (int i=0; i<number_of_particles; i++)
	{ 
		//average speed
		sum_of_velocities.x += velocity[i][0].x;
		sum_of_velocities.y += velocity[i][0].y;
	}
	velocity_average.x = sum_of_velocities.x/number_of_particles;
	velocity_average.y = sum_of_velocities.y/number_of_particles;
	for (int i=0; i<number_of_particles; i++)
	{
		velocity[i][0].x = velocity[i][0].x - velocity_average.x;
		velocity[i][0].y = velocity[i][0].y - velocity_average.y;
	}

double k=0;
for (int i=0; i<number_of_particles; i++)
{
	k+=velocity[i][0].x+velocity[i][0].y;
}
k=0;
pulse << 0 << " " << k <<endl;
	
// KIN ENERGIES 0
	for (int i=0; i<number_of_particles; i++)
	{
		kin_energy[0]+=(pow(velocity[i][0].x,2)+pow(velocity[i][0].y,2))/2;
	}	
//COORDINATES 0
		for (int i=0; i<u;  i++){
			for (int j=0; j<u; j++){
				coord[i*u+j][0].x= j*a;
				coord[i*u+j][0].y=i*h;
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
				write <<  coord[i*u+j][0].x << " " << coord[i*u+j][0].y << endl;
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
				delta[i][