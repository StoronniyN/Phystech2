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
//TRIANGLE: Lx=5 -> 0.704; Lx=7 -> 1.047    0.699 9
/*----------------------*/
string preference="triangle";
double	sigma=1.047, 
		L_x=7;
/*----------------------*/
const int n_c = 6, number_of_particles = pow(n_c,2), u=sqrt(number_of_particles); 
double L_y=L_x*(sqrt(3.0)/2.0), L=L_x*sqrt((sqrt(3.0))/2.0), a, h, V0, Vt; 
const int step = 1000; //skip 1000 coordinates
const int NTIME=1000000 /* number of steps */, control_temp=0 /* T control enable/disable 1/0 */;
int o = 0, number_of_particles_new; //T control steps
double length[number_of_particles][number_of_particles];
const double epsilon=0.0031, dt=0.001, primary_velocity = 0.001;// sigma=0.985;
double temperature_average, temperature_desired, temperature_prev, temperature, temperature_sum, kin_energy[NTIME], pot_energy[NTIME];
int stop = 0;
double MSD, MSDsum, VACF, VACFsum, dx[number_of_particles], dy[number_of_particles];
const int control_point_1=10/dt, control_point_2=50/dt, control_point_3=100/dt; //steps at which temperature control starts
double coordx[number_of_particles],
	   coordy[number_of_particles],
	   accelerationx[number_of_particles],
	   accelerationy[number_of_particles],
	   accelerationx_old[number_of_particles],
	   accelerationy_old[number_of_particles],
	   velocityx[number_of_particles],
	   velocityy[number_of_particles],
	   velocityx_old[number_of_particles],
	   velocityy_old[number_of_particles];
double random(double min, double max)
{
    return (double)(rand())/RAND_MAX*(max - min) + min;
}
double periodic_boundary_conditions_x(int i, int N, double x, double L, string preference);
double periodic_boundary_conditions_y(int i, int N, double y, double L, string preference);
double R(double d_x, double d_y, string preference, double L);
double calc_energy(double coordx[], double coordy[]);
//STRUCTURE
struct primary_pairs
{
	public:
		double x;
		double y;
};
	primary_pairs delta[number_of_particles][number_of_particles], velocity_average, sum_of_velocities;
int main(){
	if(preference=="square"){
		a=L/(n_c);
	}
	if(preference=="triangle"){
		a=L_x/(n_c);
	}
	h=a*(sqrt(3)/2);
ofstream pulse, write, kinwrite, potwrite, fullwrite, tempwrite, avtempwrite, partwrite,
MSDwrite, VACFwrite, sigma_p_energy, vacout;
write.open("data.txt");
kinwrite.open("kinenergy");
potwrite.open("potenergy");
fullwrite.open("fullenergy");
tempwrite.open("temperature");
avtempwrite.open("avtemperature");
partwrite.open("particle");
MSDwrite.open("MSD");
VACFwrite.open("VACF");
pulse.open("pulse");
sigma_p_energy.open("p_energy");
vacout.open("vacout");
	for (int i=0; i<number_of_particles; i++){
		for (int j=0; j<number_of_particles; j++){
			length[i][j]=0;
		}
	}	
	//RANDOM NUMBER GENERATOR
    srand((unsigned int)time(0));

	for (int i=0; i<number_of_particles; i++){
		for (int N=0; N<NTIME; N++){
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
		velocityx[i] =random(-primary_velocity, primary_velocity);
		//////cout << "v_x_" << i << "= " << velocity[i].x << endl;
		velocityy[i] = random(-primary_velocity, primary_velocity);
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
		sum_of_velocities.x += velocityx[i];
		sum_of_velocities.y += velocityy[i];
	}
	velocity_average.x = sum_of_velocities.x/number_of_particles;
	velocity_average.y = sum_of_velocities.y/number_of_particles;
	for (int i=0; i<number_of_particles; i++){
		velocityx[i] = velocityx[i] - velocity_average.x;
		velocityy[i] = velocityy[i] - velocity_average.y;
	}
	for (int i=0; i<number_of_particles; i++){
		velocityx_old[i]=velocityx[i];
		velocityy_old[i]=velocityy[i];
	}
double k=0;
for (int i=0; i<number_of_particles; i++)
{
	k+=velocityx[i]+velocityy[i];
}
k=0;
pulse << 0 << " " << k <<endl;
	
// KIN ENERGIES 0
	for (int i=0; i<number_of_particles; i++)
	{
		kin_energy[0]+=(pow(velocityx[i],2)+pow(velocityy[i],2))/2;
	}	
	temperature_sum+=kin_energy[0]/number_of_particles;
	temperature_average=temperature_sum;
	avtempwrite << 0 << " " << temperature_average << endl;
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
				//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
				//write <<  coordx[i*u+j] << " " << coordy[i*u+j] << endl;
			}
		}
	//cout << "" << endl;
	//cout << "" << endl;
	write << "" << endl;
	write << "" << endl;

	int n_deleted=0;
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
	
	for (int i=0; i<number_of_particles_new; i++){
			for ( int j=0; j<number_of_particles_new; j++){
				
				if (i==j)
				{
					continue;
				}
				delta[i][j].x=(coordx[i]-coordx[j]);
				delta[i][j].y=(coordy[i]-coordy[j]);
				if (preference=="square"){
					if (abs(delta[i][j].x)>L/2 )
						{
							delta[i][j].x=(delta[i][j].x/abs(delta[i][j].x))*(abs(delta[i][j].x)-L);
						}
					if (abs(delta[i][j].y)>L/2)
						{
							delta[i][j].y=(delta[i][j].y/abs(delta[i][j].y))*(abs(delta[i][j].y)-L);
						}
					length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
				}
				if (preference=="triangle"){
					if (abs(delta[i][j].x)>L_x/2 ){
							delta[i][j].x=(delta[i][j].x/abs(delta[i][j].x))*(abs(delta[i][j].x)-L_x);
						}
					if (abs(delta[i][j].y)>L_y/2){
							delta[i][j].y=(delta[i][j].y/abs(delta[i][j].y))*(abs(delta[i][j].y)-L_y);
						}
					length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
				}
			}
	}
	

	
	for (int i=0; i<number_of_particles_new; i++)
	{
		for (int j=0; j<number_of_particles_new; j++){
			if (i==j){
					continue;
			}
			accelerationx[i]+=-epsilon*(24*delta[i][j].x/pow(length[i][j],2))*( ( (pow(sigma/length[i][j],6) )) - ( (2*pow(sigma/length[i][j],12))) );
			accelerationy[i]+=-epsilon*(24*delta[i][j].y/pow(length[i][j],2))*( ( (pow(sigma/length[i][j],6) )) - ( (2*pow(sigma/length[i][j],12))) );
		}	
	}
	for (int i=0; i<number_of_particles_new; i++){
		accelerationx_old[i]=accelerationx[i];
		accelerationy_old[i]=accelerationy[i];
	}
	for (int i=0; i<number_of_particles_new; i++)
	{
		for (int j=i+1; j<number_of_particles_new; j++)
		{
			pot_energy[0]+=4*epsilon*(pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
		}
	}
  vacout << n_deleted << " "  << pot_energy[0]<< endl;



	for (int i=0; i<number_of_particles_new; i++){
		coordx[i]=coordx[i]+velocityx[i]*dt+pow(dt,2)/2*accelerationx[i];
		coordy[i]=coordy[i]+velocityy[i]*dt+pow(dt,2)/2*accelerationx[i];
		if (preference=="square"){
			if (coordx[i]>L)
			{
				coordx[i]=coordx[i]-L;
			}
			else if (coordx[i]<0)
			{
				coordx[i]=coordx[i]+L;
			}
			if (coordy[i]>L)
			{
				coordy[i]=coordy[i]-L;
			}
			else if (coordy[i]<0)
			{
				coordy[i]=coordy[i]+L;
			}
		}
		if (preference=="triangle"){

		}				
		
	}
	for (int i=0; i<number_of_particles_new; i++){
		V0=sqrt(pow(velocityx[i],2)+pow(velocityy[i],2));
	}
	kinwrite << 0 << " " << kin_energy[0] << endl;
	potwrite <<  0 << " " <<  pot_energy[0] << endl; 
	fullwrite << 0  << " " << kin_energy[0]+pot_energy[0] << endl; 

//CYCLE
for (int N=1; N<NTIME; N++){
	if (stop==1){
		cout << "ЗАВЕРШЕНИЕ" << endl;
		cout << N-1 << endl;
		break;
	}
	for (int i=0; i<number_of_particles_new; i++){
			for ( int j=0; j<number_of_particles_new; j++){
				
				if (i==j){
					continue;
				}
				delta[i][j].x=(coordx[i]-coordx[j]);
				delta[i][j].y=(coordy[i]-coordy[j]);
				
				
				//cout << "dx["<< i << "," << j<< "] = " << delta[i][j].x << endl;
				//cout << "dy["<< i << "," << j<< "] = " << delta[i][j].y << endl;
				//cout << "length["<< i << "," << j<< "] = " << length[i][j] << endl;
				if (preference=="square"){
					if (abs(delta[i][j].x)>L/2 ){
							delta[i][j].x=(delta[i][j].x/abs(delta[i][j].x))*(abs(delta[i][j].x)-L);
						}
					if (abs(delta[i][j].y)>L/2){
							delta[i][j].y=(delta[i][j].y/abs(delta[i][j].y))*(abs(delta[i][j].y)-L);
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
				accelerationx[i]+=-(24*delta[i][j].x/pow(length[i][j],2))*epsilon*( ( (pow(sigma/length[i][j],6) )) - ( (2*pow(sigma/length[i][j],12))) );
				accelerationy[i]+=-(24*delta[i][j].y/pow(length[i][j],2))*epsilon*( ( (pow(sigma/length[i][j],6) )) - ( (2*pow(sigma/length[i][j],12))) );
				//cout << "(24*delta[i][j].x/pow(length[i][j],2)) = " << (24*delta[i][j].x/pow(length[i][j],2)) << endl;
				//cout << "length[i][j] = " << length[i][j] << endl;
			}
	}
	for (int i=0; i<number_of_particles_new; i++){
		accelerationx_old[i]=accelerationx[i];
		accelerationy_old[i]=accelerationy[i];
	}
	for (int i=0; i<number_of_particles_new; i++){
		for (int j=i+1; j<number_of_particles_new; j++){
			pot_energy[N]+=4*epsilon*(pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
		}
	}
// - VELOCITIES
	for (int i=0; i<number_of_particles_new; i++){
		velocityx_old[i]=velocityx[i];
		velocityy_old[i]=velocityy[i];
		velocityx[i]=velocityx[i]+dt/2*(accelerationx_old[i] + accelerationx[i]);
		velocityy[i]=velocityy[i]+dt/2*(accelerationy_old[i] + accelerationy[i]);;
		//cout << "velocity_x_" << i << "_" << N << "= " <<  velocity[i][N].x << endl;			
		//cout << "velocity_y_" << i << "_" << N << "= " <<  velocity[i][N].y << endl;
	}

	for (int i=0; i<number_of_particles_new; i++){
		accelerationx_old[i]=0;
		accelerationx[i]=0;
		accelerationy_old[i]=0;
		accelerationy[i]=0;
	}
// KIN ENERGY AND TEMPETATURE
	for (int i=0; i<number_of_particles_new; i++){
		kin_energy[N]+=(pow(velocityx[i],2)+pow(velocityy[i],2))/2;
	}	
	temperature_sum+=kin_energy[N]/number_of_particles_new;
	temperature_average=temperature_sum/N;
	temperature=kin_energy[N]/number_of_particles_new;
	temperature_desired=0.35;
	double coeff;
	coeff=temperature_desired/temperature_average;
	if (control_temp == 1 and (N==control_point_1 or N==control_point_2 or N==control_point_3)){
		for (int i=0; i<number_of_particles_new; i++){
			velocityx[i]=sqrt(coeff)*velocityx[i];
			velocityy[i]=sqrt(coeff)*velocityy[i];
			//temperature_sum=0;
			//failback=1;
			//break;
		}
		o++;
	}
	for (int i=0; i<number_of_particles_new; i++){
			if (i==10){
				partwrite <<  coordx[i] << " " << coordy[i] << endl;
			}
			if (N % step == 0 or (N % step) % step ==0 or ((N % step) % step)%step == 0){
			write <<  coordx[i] << " " << coordy[i] << endl;
			//partwrite <<  coord[i][N].x << " " << coord[i][N].y << endl;
			if (i==10){
				//partwrite <<  coord[i][N].x << " " << coord[i][N].y << endl;
			}
		}
	}
	
	for (int i=0; i<number_of_particles_new; i++){
		coordx[i]=coordx[i]+velocityx[i]*dt+(1/2)*pow(dt,2)*accelerationx[i];
		coordy[i]=coordy[i]+velocityy[i]*dt+(1/2)*pow(dt,2)*accelerationy[i];
		if (preference=="square"){
			if (coordx[i]>L){
				coordx[i]=coordx[i]-L;
			}
			else if (coordx[i]<0){
				coordx[i]=coordx[i]+L;
			}
			if (coordy[i]>L){
				coordy[i]=coordy[i]-L;
			}
			else if (coordy[i]<0){
				coordy[i]=coordy[i]+L;
			}		
		}
		if (preference=="triangle"){
		
		}

		if ( (preference=="square" and (abs(coordx[i])>L or abs(coordy[i])>L)) or 
		     (preference=="triangle" and (abs(coordx[i])>L_x or abs(coordy[i])>L_y))
		 or kin_energy[N]-kin_energy[N-1]>1 or accelerationx[i]-accelerationx_old[i]>10 or accelerationy[i]-accelerationy_old[i]>1){
			cout << "РАСХОДИМОСТЬ" << endl;
			stop = 1;
			//cout << "velocity_x_" << i << "_" << N-1 << "= " <<  velocity[i][N-1].x << endl;			
			//cout << "velocity_y_" << i << "_" << N-1 << "= " <<  velocity[i][N-1].y << endl;
			cout << "velocity_x_" << i << "_" << N << "= " <<  velocityx[i] << endl;			
			cout << "velocity_y_" << i << "_" << N << "= " <<  velocityy[i] << endl;		
			cout << "acceleration_x_" << i << "_" << N-1 <<  " = " << accelerationx_old[i] << endl;	
			cout << "acceleration_y_" << i << "_" << N-1 <<  " = " << accelerationy_old[i] << endl;
			cout << "acceleration_x_" << i << "_" << N <<  " = " << accelerationx[i] << endl;	
			cout << "acceleration_y_" << i << "_" << N <<  " = " << accelerationy[i] << endl;
			cout << "SUM VEL = " << k << endl;
			cout << "o = " << o << endl;
			break;
		}
		dx[i]+=velocityx[i]*dt+accelerationx[i]*pow(dt,2)/2.0;
		dy[i]+=velocityy[i]*dt+accelerationy[i]*pow(dt,2)/2.0;
		//V0=sqrt(pow(velocityx_old[i],2)+pow(velocityy_old[i],2));
		MSDsum+= pow(dx[i],2)+pow(dy[i],2);/*-sqrt(pow(coord[i][0].x,2)+pow(coord[i][0].y,2))),2)*/
		Vt=sqrt(pow(velocityx[i],2)+pow(velocityy[i],2));
		VACFsum+=V0*Vt;
	}
		MSD=MSDsum/number_of_particles_new;
		VACF=VACFsum/number_of_particles_new;
		MSDwrite << N << " " << 100*MSD << endl;
		VACFwrite << N << " " << VACF << endl;
		cout << "step " << N << ": MSD = " << MSD << endl;
		cout << "step " << N << ": VACF = " << VACF << endl;
		VACFsum=0;
		MSDsum=0;
		//R=0;
		if (N % step == 0 or (N % step) % step ==0 or ((N % step) % step)%step == 0){
		write << "" << endl;
		write << "" << endl;
		}

	//total energy = kin_energyх[0]+kin_energy[0]
	velocity_average.x = 0;
	velocity_average.y = 0;
	sum_of_velocities.x = 0;
	sum_of_velocities.y = 0;
	for (int i=0; i<number_of_particles_new; i++)
	{ 
		//average speed
		sum_of_velocities.x += velocityx[i];
		sum_of_velocities.y += velocityy[i];
	}
	velocity_average.x = sum_of_velocities.x/number_of_particles_new;
	velocity_average.y = sum_of_velocities.y/number_of_particles_new;
	for (int i=0; i<number_of_particles_new; i++)
	{
		velocityx[i] = velocityx[i] - velocity_average.x;
		velocityy[i] = velocityy[i]- velocity_average.y;
	}
	for (int i=0; i<number_of_particles_new; i++)
	{
		k+=velocityx[i]+velocityy[i];
	}
	pulse << N << " " << k <<endl;
	cout << "" << endl;
	cout << "step " << N << ": kin energy = " << kin_energy[N] << endl;
	cout << "step " << N << ": pot energy = " <<  pot_energy[N] << endl;
	cout << "step " << N << ": total energy = " << kin_energy[N]+pot_energy[N] << endl;
	cout << "step " << N << ": temperature = " << temperature << endl;
	cout << "step " << N << ": temperature = " << temperature_average << endl;
	cout << "step " << N << ": pulse_sum = " << k << endl;
	k=0;
	//cout << "TOTAL energy " << N << " = " << pot_energy[N] + kin_energy[N] << endl;
	kinwrite << N << " " << setprecision(15)  <<kin_energy[N] << endl;
	potwrite <<  N << " " << setprecision(15)  <<  pot_energy[N] << endl; 
	fullwrite << N  << " " << setprecision(15) << kin_energy[N]+pot_energy[N] << endl;  
	avtempwrite << N << " " << temperature_average << endl;
	tempwrite << N << " " << temperature << endl;
	//cout << "kin energy_" << 0 << " = " << kin_energy[0]  << endl;
	//cout << "pot energy_" << 0 << " = " << pot_energy[0] << endl;
	//cout << "TOTAL energy_" << 0 << " = " << pot_energy[0] + kin_energy[0]  << endl;
	//cout << "temperature_average" << " = " << temperature_average  << endl;

	
}
write.close();
kinwrite.close();
potwrite.close();
fullwrite.close();
tempwrite.close();
avtempwrite.close();
partwrite.close();
MSDwrite.close();
VACFwrite.close();
cout << "SUM VEL = " << k << endl;
cout << "o = " << o << endl;
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

