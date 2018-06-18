#include <stdlib.h>
#include <time.h>
//#include <iostream>
//#include <limits> 
//#include <iomanip> 
#include <math.h>
//#include <fstream>
#include <stdio.h>
#include <string.h>
//#include <cstdlib>
# define M_PIl          3.141592653589793238462643383279502884L
//using namespace std;
int main(int argc, char** argv){
int rank, size;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
printf(“My rank is %d\n”, rank);
printf(“Hello, World!\n”);
MPI_Finalize();
return 0;
}
