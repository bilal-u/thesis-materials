#include <stdlib.h>
#include <omp.h>
#include <time.h>
//#inclue <ctime>
#include <stdio.h>

#include "ridc.h"

using namespace std;

class ExplicitOde : public ODE {
public:
  ExplicitOde(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt) {
    neq = my_neq;
    nt = my_nt;
    ti = my_ti;
    tf = my_tf;
    dt = my_dt;
  }
  
  void rhs(double t, double *u, double *f) {
    //  f[0]=-2.0*M_PI*sin(2.0*M_PI*t) - 2.0*( u[0]-cos(2.0*M_PI*t) );
      f[0]=u[0];
  }

  void step(double t, double * u, double * unew) {
    double* fold = new double[neq];
    rhs(t,u,fold);
  
    for (int i = 0; i < neq; i++)  {
      unew[i] = u[i] + dt*(fold[i]);
    }
    delete [] fold;
  }
};


int main(int argc, char *argv[]) {
	
 int n=400;
 
 double* tm = new double[n];
 
 for (int s=0;s<n;s++)
 {
 
 clock_t begin = clock();
 	
  int order, nt;
  double *sol;

  if (argc != 3) {
    printf("usage: <executable> <order> <nt>  >  output_file\n");
    fflush(stdout);
    exit(1);
  } else {
    order = atoi(argv[1]); // order of method
    nt = atoi(argv[2]); // number of time steps
  }

  int neq = 1;
  int ti = 0;
  int tf = 1;
  double dt = (double)(tf - ti)/nt; // compute dt
  
  // initialize ODE variable
  ExplicitOde *ode = new ExplicitOde(neq,nt,ti,tf,dt);

  sol = new double[neq];
  // specify initial condition
  
  sol[0]=1.0;
  
  // call ridc 
  ridc_fe(ode,order, sol);
 
  clock_t end = clock();
  double tim_s = double(end - begin) / CLOCKS_PER_SEC;
  
  tm[s]=tim_s;
  
  delete [] sol;
  
  }
  
  double tsum=0.0;
  
  for (int i = 0; i < n; i++)
  
      {
	  tsum=tsum+tm[i];  
	  }
  
  
  //double diff (((double)t2-(double)t1)/(double)CLOCKS_PER_SEC);
  
  printf("Time taken: %16.12fs\n",tsum/n);
  
  // output solution to screen
  // for (int i = 0; i < neq; i++){
  //  printf("%14.12f\n", sol[i]);}
   // printf("Time taken: %.12fs\n", (double)(clock() - start_s)/CLOCKS_PER_SEC); 

  
}
