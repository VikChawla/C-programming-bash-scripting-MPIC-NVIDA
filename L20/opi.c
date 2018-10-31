#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
int main(int argc, char **argv){
  double ttime = omp_get_wtime();

  long long int Ninside = 0; // number of random points inside 1/4 circle
  long long int Ntests = 1000000000;
  long long n;
  int Nthreads = atoi(argv[1]);
  omp_set_num_threads(Nthreads);
  struct drand48_data buff;

  double estpi = 0;
#pragma omp parallel num_threads(Nthreads) reduction(+ : Ninside)
  { 
    srand48_r(12345*omp_get_thread_num(), &buff);
#pragma omp for
  
  for(n=0;n<Ntests;++n){
    double x;
    double y;
    drand48_r(&buff, &x);
    drand48_r(&buff, &y);
    if(x*x+y*y<1){
      ++Ninside;
    }
  }
  }
  estpi = 4.*(Ninside/(double)Ntests);
  double time = omp_get_wtime() - stime;
  printf("estPi = %lf\n", estpi);
  printf("time=%lf\n", time);


  return 0;
}
