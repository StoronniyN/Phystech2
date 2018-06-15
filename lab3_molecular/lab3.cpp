#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <limits> 
#include <iomanip> 
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
//#include <cstdlib>
using namespace std;
/*----------------------*/
//SQUARE: Lx=5 -> 0.748; Lx=7 -> 0.985
//TRIANGLE: Lx=5 -> 0.704; Lx=7 -> 1.047
/*----------------------*/
string preference="triangle";
const double	sigma=1.047, 
				L_x=7,
				dr=0.05,
				kT=0.001;
/*----------------------*/
double dE, p;
const int n_c = 6, number_of_particles = pow(n_c,2), u=sqrt(number_of_particles); 
double a, h;
const int step = 10; //skip 1000 coordinates
const int NTIME=100000 /* number of steps */ /* T control enable/disable 1/0 */;
double length[number_of_particles][number_of_particles];
const double epsilon=0.0031;
double  L, x, y, coordx[number_of_particles],
				 coordy[number_of_particles],
				 coordx_old[number_of_particles],
				 coordy_old[number_of_particles];
int stop = 0;
double MSD, MSDsum, dx, dy,
					d_x, d_y,
					E, E_old,
					mean_E,
					mean_sq_E,
					each_E[NTIME];
double random(double min, double max){
    return (double)(rand())/RAND_MAX*(max - min) + min;
}
double periodic_boundary_conditions_x(int i, int N, double x , double L, string preference);
double periodic_boundary_conditions_y(int i, int N, double y , double L, string preference);
double R(double d_x, double d_y, string preference, double L);
double calc_energy(double coordx[], double coordy[]);
//STRUCTURE

int main(){
	if(preference=="square"){
		L=L_x*sqrt((sqrt(3.0))/2.0);
		a=L/(n_c);
	}
	if(preference=="triangle"){
		L=L_x;
		a=L/(n_c);
	}
	h=a*(sqrt(3)/2);
ofstream pulse, write, potwrite,  partwrite,
MSDwrite, VACFwrite;
write.open("data.txt");
potwrite.open("potenergy");
partwrite.open("particle");
MSDwrite.open("MSD");
VACFwrite.open("VACF");
	for (int i=0; i<number_of_particles; i++){
		for (int j=0; j<number_of_particles; j++){
			length[i][j]=0;
		}
	}	
	//RANDOM NUMBER GENERATOR
    srand((unsigned int)time(0));

	for (int i=0; i<number_of_particles; i++){
		for (int N=0; N<NTIME; N++){
			coordx[i] = 0;
			coordy[i] = 0;
		}
	}
//COORDINATES 0
		for (int i=0; i<u;  i++){
			for (int j=0; j<u; j++){
				if (preference == "triangle"){
					coordx[i*u+j]= i%2==0  ? j*a : j*a+a/2;
					coordy[i*u+j]=i*h;
				}
				if (preference == "square"){
					coordx[i*u+j]=j*a;
					coordy[i*u+j]=i*a;
				}				
				write <<  coordx[i*u+j] << " " << coordy[i*u+j] << endl;
			}
		}
	write << "" << endl;
	write << "" << endl;
	E=calc_energy(coordx, coordy);
	potwrite <<  0 << " " <<  E << endl;
	mean_E+=E;
	each_E[0]=E;
//CYCLE
for (int N=1; N<NTIME; N++){
	if (stop==1){
		cout << "ЗАВЕРШЕНИЕ" << endl;
		cout << N-1 << endl;
		break;
	}
	for (int i=0; i<number_of_particles; i++){
		for (int n=0; n<number_of_particles; n++){
			coordy_old[n]=coordy[n];
			coordx_old[n]=coordx[n];
		}
		dx=random(-dr, dr);
		coordx[i]=coordx[i]+dx;
		dy=random(-dr, dr);
		coordy[i]=coordy[i]+dy;	
		x=coordx[i];
		y=coordy[i];
		coordx[i]=periodic_boundary_conditions_x(i, N, x, L, preference);
		coordy[i]=periodic_boundary_conditions_y(i, N, y, L, preference);
		E=calc_energy(coordx, coordy);
		E_old=calc_energy(coordx_old, coordy_old);
		//cout << "E_old = " << E_old << endl;
		//cout << endl;
		dE=E-E_old;
		if (dE<=0){
			coordx[i]=coordx[i];
			MSDsum+=pow(dx,2)+pow(dy,2);
		}
		if (dE>0){
			p=random(0,1);
			if(p>exp(-(dE/kT))){
				coordx[i]=coordx_old[i];
				coordy[i]=coordy_old[i];
				E=E_old;
				}
			else{
				MSDsum+=pow(dx,2)+pow(dy,2);
			}
		}
		//cout << "MSDsum = " << MSDsum << endl;
		if (  (preference=="square" and (abs(coordx[i])>L or abs(coordy[i])>L)) or 
		     (preference=="triangle" and (abs(coordx[i])>L or abs(coordy[i])>L*(sqrt(3.0)/2.0)))){
			cout << "РАСХОДИМОСТЬ" << endl;
			stop = 1;
			cout << i << " " <<  coordx_old[i]<< " ---> " << coordx[i] << endl; 
			cout << i << " " <<  coordy_old[i] << " ---> " << coordy[i] << endl;
			break;
		}
		if (i==10){
			partwrite <<  coordx[i] << " " << coordy[i] << endl;
		}
			if (N % step == 0 or (N % step) % step ==0 or ((N % step) % step)%step == 0){
			write <<  coordx[i] << " " << coordy[i] << endl;

		}
	}
		MSD=MSDsum/number_of_particles;
		MSDsum=0;
		MSDwrite << N << " " << MSD << endl;
		cout << "step " << "MSD = " << MSD << endl;
		if (N % step == 0 or (N % step) % step ==0 or ((N % step) % step)%step == 0){
			write << "" << endl;
			write << "" << endl;
		}

	
	cout << "" << endl;
	cout << "step " << N << ": pot energy = " << E << endl;
	potwrite <<  N << " " <<  E << endl; 
	mean_E+=E;
	each_E[N]=E;
	E=0;
	E_old=0;
}
mean_E=mean_E/NTIME;
for (int N=1; N<NTIME; N++){
	mean_sq_E+=pow((each_E[N]-mean_E),2);
}
mean_sq_E=mean_sq_E/NTIME;
cout << "SIGMA E = " << mean_sq_E << endl;
write.close();
potwrite.close();
partwrite.close();
MSDwrite.close();
VACFwrite.close();

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
