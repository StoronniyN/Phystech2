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
string preference="square";
const double	sigma=0.985, 
				L_x=7,
				L=L_x*sqrt((sqrt(3.0))/2.0),
				dr=0.0001,
				kT=0.00000001;
/*----------------------*/
double dE, p, number;
const int n_c = 6, number_of_particles = pow(n_c,2), u=sqrt(number_of_particles); 
double L_y=L_x*(sqrt(3.0)/2.0), a, h, E;
const int step = 1000; //skip 1000 coordinates
const int NTIME=300000 /* number of steps */ /* T control enable/disable 1/0 */;
double energy[number_of_particles][NTIME], Epot[NTIME], E_old;
double length[number_of_particles][number_of_particles];
const double epsilon=0.0031;
double pot_energy[NTIME];
int stop = 0;
double MSD, MSDsum, VACF, VACFsum, dx, dy;
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
 	primary_pairs coord[number_of_particles][NTIME];
	primary_pairs delta[number_of_particles][number_of_particles];
int main(){
	if(preference=="square"){
		a=L/(n_c);
	}
	if(preference=="triangle"){
		a=L_x/(n_c);
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
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}
	
//COORDINATES 0
		for (int i=0; i<u;  i++){
			for (int j=0; j<u; j++){
				if (preference == "triangle"){
					coord[i*u+j][0].x= i%2==0  ? j*a : j*a+a/2;
					coord[i*u+j][0].y=i*h;
				}
				if (preference == "square"){
					coord[i*u+j][0].x=j*a;
					coord[i*u+j][0].y=i*a;
				}				
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
				write <<  coord[i*u+j][0].x << " " << coord[i*u+j][0].y << endl;
			}
		}
	//cout << "" << endl;
	//cout << "" << endl;
	write << "" << endl;
	write << "" << endl;
	for (int i=0; i<number_of_particles; i++){
			for ( int j=0; j<number_of_particles; j++){
				
				if (i==j)
				{
					continue;
				}
				delta[i][j].x=(coord[i][0].x-coord[j][0].x);
				delta[i][j].y=(coord[i][0].y-coord[j][0].y);
				if (preference=="square"){
					if (abs(delta[i][j].x)>L/2 )
						{
							delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-L);
						}
					if (abs(delta[i][j].y)>L/2)
						{
							delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-L);
						}
					length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
				}
				if (preference=="triangle"){
					if (abs(delta[i][j].x)>L_x/2 ){
							delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-L_x);
						}
					if (abs(delta[i][j].y)>L_y/2){
							delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-L_y);
						}
					length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
				}
			}
	}
	

	
	for (int i=0; i<number_of_particles; i++){
		for (int j=i+1; j<number_of_particles; j++)
		{
			energy[i][0]+=4*epsilon*(pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
			pot_energy[0]+=4*epsilon*(pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
		
		}
	}
	Epot[0]=pot_energy[0];

		//cout << "" << endl;
		//cout << "" << endl;
	potwrite <<  0 << " " <<  pot_energy[0] << endl;
	
	for (int i=0; i<number_of_particles; i++){
		dx=random(-dr, dr);
		coord[i][1].x=coord[i][0].x+dx;
		dy=random(-dr, dr);
		coord[i][1].y=coord[i][0].y+dy;
		
	}
E_old=pot_energy[0];
//CYCLE
for (int N=1; N<NTIME; N++){
	if (stop==1){
		cout << "ЗАВЕРШЕНИЕ" << endl;
		cout << N-1 << endl;
		break;
	}
	for (int i=0; i<number_of_particles; i++){
			for ( int j=0; j<number_of_particles; j++){
				
				if (i==j){
					continue;
				}
				delta[i][j].x=(coord[i][N].x-coord[j][N].x);
				delta[i][j].y=(coord[i][N].y-coord[j][N].y);
				if (preference=="square"){
					if (abs(delta[i][j].x)>L/2 ){
							delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-L);
						}
					if (abs(delta[i][j].y)>L/2){
							delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-L);
						}
				}
				if (preference=="triangle"){
					if (abs(delta[i][j].x)>L_x/2 ){
							delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-L_x);
						}
					if (abs(delta[i][j].y)>L_y/2){
							delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-L_y);
						}
				}
				length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
			}
		for (int j=i+1; j<number_of_particles; j++){
			energy[i][N]+=4*epsilon*(pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
			//pot_energy[N]+=4*epsilon*(pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
		}
		Epot[N]=0;
		for (int q=0; q<number_of_particles; q++){
			for (int w=q+1; w<number_of_particles; w++){
				Epot[N]+=4*epsilon*(pow((sigma/length[q][w]),12) - pow((sigma/length[q][w]),6) );
			}
		}
		
		
		//dE=energy[i][N]-energy[i][N-1];
		//dE=pot_energy[N]-pot_energy[N-1];
		dE=Epot[N]-E_old;
		if (dE<=0){
			coord[i][N]=coord[i][N];
			MSDsum+=pow(dx,2)+pow(dy,2);
			E+=dE;
		}
		if (dE>0){
			p=random(0,1);
			if(p>exp(-(dE/kT))){
				coord[i][N].x=coord[i][N-1].x;
				coord[i][N].y=coord[i][N-1].y;
				pot_energy[N]=pot_energy[N-1];
				//pot_energy[N]-=dE;
				
				}
			else{
				MSDsum+=pow(dx,2)+pow(dy,2);
				E+=dE;
			}
		E_old=Epot[N];
		}
		//pot_energy[N]+=energy[i][N];
		//cout << "dE = " << dE << endl;
		//cout << "Epot["<< i << "] = " << Epot[i] << endl;
		//energy[i][N]=0;
		//pot_energy[N]=energy[i][N];
	}
	pot_energy[N]=E;

	
	for (int i=0; i<number_of_particles; i++){
		if (  (preference=="square" and (abs(coord[i][N+1].x)>L or abs(coord[i][N+1].y)>L)) or 
		     (preference=="triangle" and (abs(coord[i][N+1].x)>L_x or abs(coord[i][N+1].y)>L_y))){
			cout << "РАСХОДИМОСТЬ" << endl;
			stop = 1;
			cout << i << " " <<  coord[i][N-1].x << " ---> " << coord[i][N].x << endl; 
			cout << i << " " <<  coord[i][N-1].y << " ---> " << coord[i][N].y << endl;
			for (int i=0; i<u;  i++){
				for (int j=0; j<u; j++){
					//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
					cout <<  coord[i*u+j][N].x << " " << coord[i*u+j][N].y << endl;
				}
			}
			break;
		}
		if (i==10){
			partwrite <<  coord[i][N].x << " " << coord[i][N].y << endl;
		}
			//write <<  coord[i][N].x << " " << coord[i][N].y << endl;
			
			if (N % step == 0 or (N % step) % step ==0 or ((N % step) % step)%step == 0){
			write <<  coord[i][N].x << " " << coord[i][N].y << endl;

		}
		//cout << coord[i][N].x << endl;
		//cout << coord[i][N].y << endl;
		//R+= pow((coord[i][N].x-coord[i][N-1].x),2)+pow((coord[i][N].y-coord[i][N-1].y),2);
		//dx[i]+=velocity[i][N].x*dt+acceleration[i][N].x*pow(dt,2)/2.0;
		//dy[i]+=velocity[i][N].y*dt+acceleration[i][N].y*pow(dt,2)/2.0;
		/*-sqrt(pow(coord[i][0].x,2)+pow(coord[i][0].y,2))),2)*/
		//Zsum+=V0*sqrt(pow(velocity[i][N].x,2)+pow(velocity[i][N].y,2));
	}
		MSD=MSDsum/number_of_particles;
		VACF=VACFsum/number_of_particles;
		MSDwrite << N << " " << MSD << endl;
		cout << "step " << "MSD = " << MSD << endl;
		cout << "step " << "VACF = " << VACF << endl;
		//R=0;
		if (N % step == 0 or (N % step) % step ==0 or ((N % step) % step)%step == 0){
		write << "" << endl;
		write << "" << endl;
		}

	
	cout << "" << endl;
	cout << "step " << N << ": pot energy = " << pot_energy[N] << endl;
	//cout << "TOTAL energy " << N << " = " << pot_energy[N] + kin_energy[N] << endl;
	potwrite <<  N << " " <<  pot_energy[N] << endl; 
	//cout << "kin energy_" << 0 << " = " << kin_energy[0]  << endl;
	//cout << "pot energy_" << 0 << " = " << pot_energy[0] << endl;
	//cout << "TOTAL energy_" << 0 << " = " << pot_energy[0] + kin_energy[0]  << endl;
	//cout << "temperature_average" << " = " << temperature_average  << endl;

	for (int i=0; i<number_of_particles; i++){
		dx=random(-dr, dr);
		coord[i][N+1].x=coord[i][N].x+dx;
		dy=random(-dr, dr);
		coord[i][N+1].y=coord[i][N].y+dy;
		if (preference=="square"){
			if (coord[i][N+1].x>L)
			{
				coord[i][N+1].x=0;
			}
			else if (coord[i][N+1].x<0)
			{
				coord[i][N+1].x=L;
			}
			if (coord[i][N+1].y>L)
			{
				coord[i][N+1].y=0;
			}
			else if (coord[i][N+1].y<0)
			{
				coord[i][N+1].y=L;
			}		
		}
		if (preference=="triangle"){
			if (coord[i][N+1].x>L_x)
			{
				coord[i][N+1].x=coord[i][N+1].x-L_x;
			}
			else if (coord[i][N+1].x<0)
			{
				coord[i][N+1].x=coord[i][N+1].x+L_x;
			}
			if (coord[i][N+1].y>L_y)
			{
				coord[i][N+1].y=coord[i][N+1].y-L_y;
			}
			else if (coord[i][N+1].y<0)
			{
				coord[i][N+1].y=coord[i][N+1].y+L_y;
			}		
		}
	}
}
write.close();
potwrite.close();
partwrite.close();
MSDwrite.close();
VACFwrite.close();

	return 0;
}
