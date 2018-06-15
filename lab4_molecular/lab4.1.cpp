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
# define M_PIl          3.141592653589793238462643383279502884L
using namespace std;
/*----------------------*/
//SQUARE: Lx=5 -> 0.748; Lx=7 -> 0.985
//TRIANGLE: Lx=5 -> 0.704; Lx=7 -> 1.047
/*----------------------*/
string preference="triangle";
string magnetic_preference="FM";
double 			J0=0.01;
const double	sigma=1.047, 
				L_x=7,
				dr=0.001,
				kT=0.0000001,
				r_c=2,
				d_phi=0.001;
/*----------------------*/
double dE, p;
const int n_c = 6, number_of_particles = pow(n_c,2), u=sqrt(number_of_particles); 
double a, h;
const int step = 100; //skip 1000 coordinates
const int NTIME=100000 /* number of steps */ /* T control enable/disable 1/0 */;
double length[number_of_particles][number_of_particles];
const double epsilon=0.0031;
double  L, x, y, coordx[number_of_particles],
				 coordy[number_of_particles],
				 coordx_old[number_of_particles],
				 coordy_old[number_of_particles],
				 phi[number_of_particles],
				 phi_old[number_of_particles],
				 theta[number_of_particles],
				 theta_old[number_of_particles];

int stop = 0, number_of_particles_new, deleted_particles=5;
double MSD, MSDsum, dx, dy,
					d_x, d_y,
					E, E_old;
double random(double min, double max){
    return (double)(rand())/RAND_MAX*(max - min) + min;
}
int random_int(int min, int max){
    return (int)(rand())/RAND_MAX*(max - min) + min;
}
double periodic_boundary_conditions_x(int i, int N, double x , double L, string preference);
double periodic_boundary_conditions_y(int i, int N, double y , double L, string preference);
double R(double d_x, double d_y, string preference, double L);
double calc_energy(double coordx[], double coordy[], double phi[], double theta[]);
//STRUCTURE

int main(){
	if(magnetic_preference=="AFM"){
		J0=-J0;
	}
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
MSDwrite, VACFwrite, spinwrite;
write.open("data.txt");
potwrite.open("potenergy");
partwrite.open("particle");
MSDwrite.open("MSD");
VACFwrite.open("VACF");
spinwrite.open("spin");
	for (int i=0; i<number_of_particles; i++){
		for (int j=0; j<number_of_particles; j++){
			length[i][j]=0;
		}
	}	
	//RANDOM NUMBER GENERATOR
    srand((unsigned int)time(0));

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
				
			}
		}
		for (int i=0; i<number_of_particles;  i++){
			phi[i]=random(-M_PI, +M_PI);
			theta[i]=random(-M_PI, +M_PI);
		}		
	int n_deleted=0;//rand()%(number_of_particles/4)+1;
	//number_of_particles_new=number_of_particles-n_deleted;
	int i_deleted[n_deleted];
	for (int i=0; i<n_deleted;  i++){
		i_deleted[i]=rand()%(number_of_particles)+1;
	}
	
	for (int i=0; i<n_deleted;  i++){
		for (int j=i+1; j<n_deleted;  j++){
			if (i_deleted[i]==i_deleted[j]){
				int tmp=i_deleted[j];
				i_deleted[j]=i_deleted[n_deleted];
				i_deleted[n_deleted]=tmp;
				n_deleted=n_deleted-1;
				
			}	
		}
	}	
	number_of_particles_new=number_of_particles-n_deleted;
	
	for (int i=0; i<n_deleted;  i++){
		for (int j=i+1; j<n_deleted;  j++){
			if (i_deleted[j]<i_deleted[i]){
				int tmp=i_deleted[i];
				i_deleted[i]=i_deleted[j];
				i_deleted[j]=tmp;
			}
		}
	}
	for (int j=0; j<n_deleted; j++){
		for (int i=0; i<number_of_particles;  i++){
			if (i>=i_deleted[j]){
				coordx[i]=coordx[i+1];
				coordy[i]=coordy[i+1];
				//i_deleted[j+1]=i_deleted[j+1]+1;
			}
		}
	}
	for (int i=0; i<number_of_particles_new;  i++){
		write <<  coordx[i] << " " << coordy[i] << endl;
		spinwrite <<  coordx[i] << " " << coordy[i]<<  " " << sin(theta[i])*sin(phi[i]) << " " << sin(theta[i])*cos(phi[i]) << endl;
	}
	
	write << "" << endl;
	write << "" << endl;
	E=calc_energy(coordx, coordy, phi, theta);
	potwrite <<  0 << " " <<  E << endl;
//CYCLE
for (int N=1; N<NTIME; N++){
	if (stop==1){
		cout << "ЗАВЕРШЕНИЕ" << endl;
		cout << N-1 << endl;
		break;
	}
	for (int i=0; i<number_of_particles_new; i++){

		for (int n=0; n<number_of_particles_new; n++){
			coordy_old[n]=coordy[n];
			coordx_old[n]=coordx[n];
			phi_old[n]=phi[n];
			theta_old[n]=theta[n];
		}
		dx=random(-dr, dr);
		coordx[i]=coordx[i]+dx;
		dy=random(-dr, dr);
		coordy[i]=coordy[i]+dy;
		phi[i]+=random(-d_phi*M_PI,d_phi*M_PI);
		theta[i]+=random(-d_phi*M_PI,d_phi*M_PI);
		x=coordx[i];
		y=coordy[i];
		coordx[i]=periodic_boundary_conditions_x(i, N, x, L, preference);
		coordy[i]=periodic_boundary_conditions_y(i, N, y, L, preference);
		E=calc_energy(coordx, coordy, phi, theta);
		E_old=calc_energy(coordx_old, coordy_old, phi_old, theta_old);
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
				phi[i]=phi_old[i];
				theta[i]=theta_old[i];
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
			spinwrite <<  coordx[i] << " " << coordy[i] <<  " " << sin(theta[i])*sin(phi[i]) << " " << sin(theta[i])*cos(phi[i]) << endl;

		}
	}
		MSD=MSDsum/number_of_particles_new;
		MSDsum=0;
		MSDwrite << N << " " << MSD << endl;
		cout << "step " << "MSD = " << MSD << endl;
		if (N % step == 0 or (N % step) % step ==0 or ((N % step) % step)%step == 0){
			write << "" << endl;
			write << "" << endl;
			spinwrite << "" << endl;
			spinwrite << "" << endl;
		}

	
	cout << "" << endl;
	cout << "step " << N << ": pot energy = " << E << endl;
	potwrite <<  N << " " <<  E << endl; 
	E=0;
	E_old=0;
}
write.close();
potwrite.close();
partwrite.close();
MSDwrite.close();
VACFwrite.close();
	for (int i=0; i<n_deleted;  i++){
		cout << i_deleted[i] << endl;
	}	
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
double calc_energy(double coordx[], double coordy[], double phi[], double theta[]){
	double energy=0;
	double hevis;
	for (int i=0; i<number_of_particles_new; i++){
			for ( int j=0; j<number_of_particles_new; j++){
				if (i==j){
					continue;
				}
				d_x=(coordx[i]-coordx[j]);
				d_y=(coordy[i]-coordy[j]);
				length[i][j]=R(d_x, d_y, preference, L);
			}
		for (int j=i+1; j<number_of_particles_new; j++){
			hevis= r_c - length[i][j]<0 ? 0 : 1;
					energy+=4*epsilon*(pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) )-J0*pow((1-length[i][j]/r_c),3)*hevis*(sin(theta[i])*sin(theta[j])*sin(phi[i])*sin(phi[j])+sin(theta[i])*sin(theta[j])*cos(phi[i])*cos(phi[j])+cos(theta[j]));
		}
	}
	return(energy);
}

