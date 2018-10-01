#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv){

  long long int Ninside = 0; // number of random points inside 1/4 circle
  long long int Ntests = 1000000000;
  long long n;

 int rank, size;
  MPI_Init(&argc, &argv);
  int messageN = size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &size);

 MPI_Status status;
int *messageIn = (int*) malloc(messageN*sizeof(int));
int messageSource = rank - 1;
int messageDest = rank+1;
int messageTag = 999;
 long long int sum = 0;
  
  double estpi = 0;

  srand48(12345);


  
  for(n=0;n<Ntests;++n){
    double x = drand48();
    double y = drand48();
    
    if(x*x+y*y<1)
      {
      ++Ninside;
    }
  }

  //if(rank == 0)
  //{
  //	for(int i = 1; i < size; i++)
  //	{
		
  //	MPI_Recv(messageIn, 1, MPI_LONG_LONG_INT, i,  messageTag,MPI_COMM_WORLD, &status);
		
  //	}
  // }

  //if(rank != 0)
  //{
	
	
  //	MPI_Send(&Ninside, 1, MPI_LONG_LONG_INT, 0, messageTag,	MPI_COMM_WORLD);
		

  //}

  MPI_Reduce(&Ninside, &sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(rank == 0)
    {
      estpi = 4.*(Ninside/(double)Ntests);
      printf("estPi = %lf\n", estpi);
    }
  
MPI_Finalize();
  return 0;


}
