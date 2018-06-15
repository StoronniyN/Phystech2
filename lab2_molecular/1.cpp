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
const double sigma=1, epsilon=1;
double U, r;
int main (){
	ofstream writelg;
	writelg.open("lgpotential");
	r=0.00001;
	while(r<10000){
	U=4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));
	r+=0.001;
	writelg << r << " " << U << endl;
	}
	return 0;
}

