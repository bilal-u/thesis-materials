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
     
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Mass = mass matrix /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void mass(double t, double *y, double **M)
   
   {
     
    M[0][0]=4.0;
    M[0][1]=-1.0;
    M[1][0]=-1.0;
    M[1][1]=4.0;
         
   }  
   
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// gauss elimination //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void gauss(double *F, double **J, double *xn)

  {
   
   // creating augmented matrix

    vector< vector<double> > A(neq,vector<double>(neq+1));
 
      // Read input data

             for (int j=0; j<neq; j++)
                {
              for (int k=0; k<neq; k++)
                  {
                 A[j][k]=J[j][k];
                  }
                }

            for (int j=0; j<neq; j++)
                {
                A[j][neq]=F[j];
                }

      // solving AX=B 

   for (int i=0; i<neq; i++)
      {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<neq; k++)
          {
            if (abs(A[k][i]) > maxEl)
            {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
          }

        // Swap maximum row with current row (column by column)

        for (int k=i; k<neq+1;k++)
           {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
           }

        // Make all rows below this one 0 in current column

        for (int k=i+1; k<neq; k++)
        {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<neq+1; j++)
            {
                if (i==j)
                {
                    A[k][j] = 0;
                } else
                   {
                    A[k][j] += c * A[i][j];
                   }
            }
        }
     }  // i loop

    // Solve equation Ax=b for an upper triangular matrix A

    vector<double> x(neq);

   
    for (int i=neq-1; i>=0; i--)
      {
        x[i] = A[i][neq]/A[i][i];
        for (int k=i-1;k>=0; k--)
        {
            A[k][neq] -= A[k][i] * x[i];
        }
      }

    for (int k=0;k<neq; k++)
        {
            xn[k] = x[k];
        }

}      

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////// rhs of the ode y'=f(t,y) ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

   void rhs(double t, double *u, double *f)
   
    {
			
	double *g = new double[neq];					

        g[0]=u[0] + 4*u[1];
        g[1]=-4*u[0] - u[1]; 
        
        double **L = new double*[neq];

         for (int j=0;j<neq;j++) 
             {
             L[j] = new double[neq];
             }
        
        mass(t,u,L);
        
        double *w = new double[neq];

     // linear solve

      gauss(g,L,w);

    
       for (int i=0;i<neq;i++)

         {
         f[i]=w[i];
         }

     for (int i=0; i<neq; i++)
      {
       delete [] L[i];
      }
  
    delete [] L;

    delete [] g;
    
    delete [] w;
    
    
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

void jac(double t, double *x, double *xold, double **J)
  {

  double dx=0.00000001;

   
  for (int i=0;i<neq;i++)

     {

     vector<double> xx(neq);

      for (int j=0;j<neq;j++)

       {
       xx[j]=x[j];
       }

     xx[i]=x[i]+dx; 

    double *xn = new double[neq];
    
    for (int j=0;j<neq;j++)

       {
       xn[j]=xx[j];
       }
     
    double *fx = new double[neq];

    fun(t,x,xold,fx);

   double *fxx = new double[neq];

    fun(t,xn,xold,fxx);

     for (int k=0;k<neq;k++)
      
        {
         J[k][i]=(fxx[k]-fx[k])/dx;
        }

    delete [] fx;
   
    delete [] xn;

    delete [] fxx;

    }
 

 }

////////////////////////////////////////////////////////////////////////////////
////////////////////////// l2 norm of vector ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    void l2_norm(double *w, int n,double *norm)
        {
         double accum = 0.;
         for (int i = 0; i < n; ++i)
          {
            accum += w[i] * w[i];
          }
         norm[0]=sqrt(accum);

         }

////////////////////////////////////////////////////////////////////////////////
////////////////////// diff = difference of two vectors ////////////////////////
////////////////////////////////////////////////////////////////////////////////

   void diff(double *yn, double *xn, double *vn)
           {
           for (int i=0;i<neq;i++)
            {
            vn[i]=(yn[i]-xn[i]);
            }
           }
      
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Newton's Solver ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//void newton(int max_iter, double tol, double t, double *u, double *unew)
void step(double t, double *u, double *unew)

{   
       //  printf("step being called \n");
      
         int max_iter=100;

        double tol=0.00000001;
	  int it_count=0;
	  vector<double> p(neq);
      
      vector<double> z(neq);
        
          for (int j=0;j<neq;j++)
              {
               p[j]=u[j];
              }
              
        double *w = new double[neq];   
        
        for (int j=0;j<neq;j++)
              {
               w[j]=u[j];
              }
     // Newton's iteration starts here!

          for (int l=0;l<max_iter;l++)
            {
		    if (l==max_iter)
                {
			    printf("not converged\n");
			    } 		
				
		    it_count=it_count+1;		
			     
             double *v = new double[neq];

             for (int j=0;j<neq;j++)
                 {
                 v[j]=p[j];
                 }
             
             double *Fn = new double[neq];

              // calling function F(X)=0
            
             fun(t,v,w,Fn);
             
             double **Jn = new double*[neq];

               for (int j=0;j<neq;j++) 
                  {
                  Jn[j] = new double[neq];
                  }
             
              // calling jacobian

              jac(t,v,w,Jn); 

         double *x=new double[neq];

         // calling linear solver for x=Jn/Fn

         gauss(Fn,Jn,x); 
         
         for (int j=0;j<neq;j++)

             {
              z[j]=v[j]-x[j];
             }
             
        double *dv=new double[neq];

         for (int j=0;j<neq;j++)

            {
            dv[j]=abs(z[j]-v[j]);
            }

       double *norm=new double[1];
       
       l2_norm(dv,neq,norm);
       
              
          if (norm[0]<tol)

            {
		  //  printf("converged\n");		
		    break;
            }
        
            p=z; 

         delete [] x;

         delete [] Fn;

          delete [] v;
          
          delete [] dv;
          
          delete [] norm;

          for (int i=0; i<neq; i++)
            {
            delete [] Jn[i];
            }
           delete [] Jn;
      } // end newton loop l
       
            
       for (int j=0;j<neq;j++)
         {
          unew[j]=z[j];
         }      
      
       
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

  int neq = 2;
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
  
    sol[0]=0.0;

    sol[1]=1.0;
  

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


