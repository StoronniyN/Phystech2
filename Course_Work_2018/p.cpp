#include "mpi.h"
#include <iostream>
//#include <windows.h>
using namespace std;

double Fact(int n)
{
  if (n==0)
   return 1;
  else 
   return n*Fact(n-1);
}

int main(int argc, char *argv[])
{
  //SetConsoleOutputCP(1251);
  int n;
  int myid;
  int numprocs;
  int i;
  int rc;
  long double drob,drobSum=0,Result, sum;
  double startwtime = 0.0;
  double endwtime;
 
  n =argv[1];

  if (rc= MPI_Init(&argc, &argv)) 
  { 
    cout << "Ошибка запуска, выполнение остановлено " << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  } 

MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
MPI_Comm_rank(MPI_COMM_WORLD,&myid); 

 
  if (myid == 0)
  { 
  
    startwtime = MPI_Wtime();
  }
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  for (i = myid; i <= n; i += numprocs)
  { 
    drob = 1/Fact(i);
    drobSum += drob;
  }

MPI_Reduce(&drobSum, &Result, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
cout.precision(20);
  if (myid == 0)
  {   
    cout << Result << endl; 
    endwtime = MPI_Wtime();
    cout << (endwtime-startwtime)*1000 << endl;     
  }

MPI_Finalize();
return 0;
}
