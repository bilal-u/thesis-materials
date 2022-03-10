#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <vector>
#include <stdio.h>
#include <math.h>

#include "ridc.h"

using namespace std;

class ImplicitOde : public ODE
   {
  public:
  ImplicitOde(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt) 
     {
    neq = my_neq;
    nt = my_nt;
    ti = my_ti;
    tf = my_tf;
    dt = my_dt;
     }
     
/////////////////////////////////////////////////////////////////////////////////////
//////////////////////// rhs of the ode y'=f(t,y) ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

   void rhs(double t, double *u, double *f)
   
    {
    f[0]=-2.0*M_PI*sin(2.0*M_PI*t) - 2.0*( u[0]-cos(2.0*M_PI*t) );
   // f[0]=u[0];
    }

/////////////////////////////////////////////////////////////////////////////////////
////////////////////  function in the form F(x)=0  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    void fun(double t, double *p, double *u, double *Fn)

      {

        double* frh = new double[neq];

           rhs(t,p,frh);
   
        for (int j=0;j<neq;j++)

               {
               Fn[j]=p[j]-u[j]-dt*frh[j];
               }

      delete [] frh; 

     } 

////////////////////////////////////////////////////////////////////////////////
//////////////////// Numerical jacobian of F(x)=0 //////////////////////////////
////////////////////////////////////////////////////////////////////////////////



//  J[0]=1.0-2*dt;
  
void jac(double t, double *x, double *xold, double *J)
  {

  J[0]=1.0-dt*( -2.0*M_PI*sin(2.0*M_PI*t) - 2.0*( 1.0-cos(2.0*M_PI*t) ) );
 

  }

////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Newton's Solver ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void step(double t, double *u, double *unew)

{   
	   int max_iter=100;

        double tol=0.00000001;
        
	   int it_count=0;
	  
	   double p,z;
       
        double *w = new double[neq];
        
        w[0]=u[0];
       
        p=u[0];
        
       // Newton's iteration starts here!

          for (int l=0;l<max_iter;l++)
            {
		    			
		    it_count=it_count+1;		
			     
             double *v = new double[neq];

             v[0]=p;
             
             double *Fn = new double[neq];

              // calling function F(X)=0
            
             fun(t,v,w,Fn);
             
             double *Jn = new double[neq];

              // calling jacobian

              jac(t,v,w,Jn); 

              z=v[0]-Fn[0]/Jn[0];
     
             if (abs(z-v[0])<tol)

            {
		    break;
            }
        
            p=z; 

        // delete [] x;

         delete [] Fn;

          delete [] v;
      
      } // end newton loop l
       
            
 unew[0]=z;
 
 delete [] w;
      
 } // newton
 
};


int main(int argc, char *argv[]) 

 {
  int start_s=clock(); 
  int order, nt;
  double *sol;

  if (argc != 3)
     {
    printf("usage: <executable> <order> <nt>  >  output_file\n");
    fflush(stdout);
    exit(1);
    } 
   else 
      {
      order = atoi(argv[1]); // order of method
      nt = atoi(argv[2]); // number of time steps
       }

  int neq = 1;
  int ti = 0;
  int tf = 1;
/////////////////// for order of accuracy ///////////////////////
/////////////////// number of step is double //////////////////////
////// l2 norm of the errors is stored in er_st vecttor////////////
///////////////////////////////////////////////////////////////////

  
  double dt = (double)(tf - ti)/double(nt); // compute dt

  
  // initialize ODE variable

  ImplicitOde *ode = new ImplicitOde(neq,nt,ti,tf,dt);

  sol = new double[neq];

  // specify initial condition
  
    sol[0]=1.0;

 // call ridc 

  ridc_be(ode, order, sol);

  // output solution to screen

  for (int i = 0; i < neq; i++)
     {
    printf("%16.14f\n", sol[i]);
     }
delete [] sol;
printf("Time taken: %.12fs\n", (double)(clock() - start_s)/CLOCKS_PER_SEC); 

}


