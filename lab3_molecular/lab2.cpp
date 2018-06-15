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
const int n_c = 9, number_of_particles = pow(n_c,2), u=sqrt(number_of_particles); 
double  a ,h, L, L_x=7, L_y=sqrt(3)*L_x/2, coordx[number_of_particles], coordy[number_of_particles];
const int step = 1000;//skip 1000 coordinates
const int NTIME=1000000 /* number of steps */, control_temp=0 /* T control enable/disable 1/0 */;
int o = 0; //T control steps
double length[number_of_particles][number_of_particles];
const double epsilon=0.0031, lattice_constant=40, dt=0.001, primary_velocity = 1;
double temperature_average, temperature_desired, temperature_prev, temperature, temperature_sum, kin_energy[NTIME], pot_energy[NTIME];
int stop = 0;
double  dx, dy, sigma=0.01, d_x, d_y, x, y;
const int control_point_1=10/dt, control_point_2=50/dt, control_point_3=100/dt; //steps at which temperature control starts
double random(double min, double max)
{
    return (double)(rand())/RAND_MAX*(max - min) + min;
}
double calc_energy(double coordx[], double coordy[]);
double R(double d_x, double d_y, string preference, double L);
double periodic_boundary_conditions_x(int i, int N, double x , double L, string preference);
double periodic_boundary_conditions_y(int i, int N, double y , double L, string preference);
//STRUCTURE
struct primary_pairs
{
	public:
		double x;
		double y;
};
	primary_pairs delta[number_of_particles][number_of_particles];
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
			coordx[i] = 0;
			coordy[i] = 0;
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
			coordx[i] = 0;
			coordy[i] = 0;
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
		L=L_x*sqrt((sqrt(3.0))/2.0);
		for (int i=0; i<u;  i++)
		{
			for (int j=0; j<u; j++){
				coordx[i*u+j]= j*a;
				coordy[i*u+j]=i*a;
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
				if(preference=="square5"){
				write_sq_5 <<  coordx[i*u+j] << " " << coordy[i*u+j] << endl;
				}
				if(preference=="square7"){
				write_sq_7 <<  coordx[i*u+j] << " " << coordy[i*u+j] << endl;
				}
			}
		}

	}
		//cout << "" << endl;
		//cout << "" << endl;
		write << "" << endl;
		write << "" << endl;
		while (sigma < 2){
			pot_energy[0]+=calc_energy(coordx, coordy);
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
			coordx[i] = 0;
			coordy[i] = 0;
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
			coordx[i] = 0;
			coordy[i] = 0;
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
				coordx[i*u+j]= j*a;
				coordy[i*u+j]=i*a;
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
				if(preference=="square7"){
				write_sq_7 <<  coordx[i*u+j] << " " << coordy[i*u+j] << endl;
			}
			}
		}

	}
		//cout << "" << endl;
		//cout << "" << endl;
		write << "" << endl;
		write << "" << endl;
		
		while (sigma < 2){
			pot_energy[0]+=calc_energy(coordx, coordy);
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
			coordx[i] = 0;
			coordy[i] = 0;
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
			coordx[i] = 0;
			coordy[i] = 0;
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
				coordx[i*u+j]= i%2==0  ? j*a : j*a+a/2;
				coordy[i*u+j]=i*h;
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
				if(preference=="triangle5"){
				write_tr_5 <<  coordx[i*u+j] << " " << coordy[i*u+j] << endl;
			}
			}
		}

	}
		//cout << "" << endl;
		//cout << "" << endl;
		write << "" << endl;
		write << "" << endl;

		while (sigma < 2){
			pot_energy[0]+=calc_energy(coordx, coordy);
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
			coordx[i] = 0;
			coordy[i] = 0;

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
			coordx[i] = 0;
			coordy[i] = 0;
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
				coordx[i*u+j]= i%2==0  ? j*a : j*a+a/2;
				coordy[i*u+j]=i*h;
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
			write_tr_7 <<  coordx[i*u+j] << " " << coordy[i*u+j] << endl;
			}
		}

	}
		//cout << "" << endl;
		//cout << "" << endl;
		write << "" << endl;
		write << "" << endl;
		while (sigma < 2){
			pot_energy[0]+=calc_energy(coordx, coordy);

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
double periodic_boundary_conditions_x(int i, int N, double x, double L, string preference){
	if (preference=="triangle"){
		if (x>L){
			x=0;
		}
		else if (x<0){
			x=L;
		}
	}
	if (preference=="square"){
		if (x>L){
			x=0;
		}
		else if (x<0){
			x=L;
		}
	}
	return (x);
}
double periodic_boundary_conditions_y(int i, int N, double y, double L, string preference){
	if (preference=="triangle"){
		if (y>L*(sqrt(3.0)/2.0)){
			y=0;
		}
		else if (y<0){
			y=L*(sqrt(3.0)/2.0);
		}
	}
	if (preference=="square"){
		if (y>L){
			y=0;
		}
		else if (y<0){
			y=L;
		}
	}
	return (y);
}
double R(double d_x, double d_y, string preference, double L){
	double length;
	if (preference=="square"){
		if (abs(d_x)>L/2){
			d_x=((d_x)/(abs(d_x)))*(abs(d_x)-L);
		}
		if (abs(d_y)>L/2){
			d_y=((d_y)/(abs(d_y)))*(abs(d_y)-L);
		}
	}
	if (preference=="triangle"){
		if (abs(d_x)>L/2 ){
			d_x=((d_x)/(abs(d_x)))*(abs(d_x)-L);
		}
		if (abs(d_y)>L*(sqrt(3.0)/2.0)/2){
			d_y=((d_y)/(abs(d_y)))*(abs(d_y)-L*(sqrt(3.0)/2.0));
		}
	}
	length=sqrt(pow(d_x,2)+pow(d_y,2));
	return(length);
}

double calc_energy(double coordx[], double coordy[]){
	double energy=0;
	for (int i=0; i<number_of_particles; i++){
			for ( int j=0; j<number_of_particles; j++){
				if (i==j){
					continue;
				}
				d_x=(coordx[i]-coordx[j]);
				d_y=(coordy[i]-coordy[j]);
				length[i][j]=R(d_x, d_y, preference, L);
			}
		for (int j=i+1; j<number_of_particles; j++){
			energy+=4*epsilon*(pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
		}
	}
	return(energy);
}
