#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <vector>
#include <stdio.h>
#include <math.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_vector.h>


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
////////////////////// diff = difference of two vectors ////////////////////////
////////////////////////////////////////////////////////////////////////////////

void diff(double *yn, double *xn, double *vn,int N)
   {
     for (int i=0;i<N;i++) {
        vn[i]=(yn[i]-xn[i]); }
    }

////////////////////////////////////////////////////////////////////////////////
////////////////////// Rho = Mesh density function /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void Rho(double *y, double *rho_sm)
 
  {
   int N=neq/2;
        
   vector<double> u(N);
   vector<double> x(N);
      
   for (int j=0;j<N;j++) {
       u[j] =y[j];
       x[j] =y[j+N]; } 

   vector<double> rho(N);

   vector<double> v(N);

   for (int j=1;j<N-1;j++) {
      v[j]=2.0/(x[j+1]-x[j-1])*( (u[j+1]-u[j])/(x[j+1]-x[j])-(u[j]-u[j-1])/(x[j]-x[j-1]) ); }
      
      v[0] = 2.0*((x[1]-x[0])*(u[2]-u[0])-(x[2]-x[0])*(u[1]-u[0]))/((x[2]-x[0])*(x[1]-x[0])*(x[2]-x[1]));
      v[N-1] = 2.0*((x[N-2]-x[N-1])*(u[N-3]-u[N-1])-(x[N-3]-x[N-1])*(u[N-2]-u[N-1]))/((x[N-3]-x[N-1])*(x[N-2]-x[N-1])*(x[N-3]-x[N-2]));
  
   for (int j=0;j<N;j++) {
	   rho[j]=rho[j]+pow(v[j],2); }
        
   // alpha calculation
   
   double gamma = 1.0/3.0;
   
   double Alpha =0.0;  

   for (int j=1;j<N;j++) {
      Alpha=Alpha+ 0.5*( pow(rho[j],gamma) + pow(rho[j-1],gamma) )*(x[j]-x[j-1]); }

    Alpha = pow(Alpha,3);

    vector<double> rh(N);

    for (int j=0;j<N;j++) {
       rh[j]=pow( (1.0+(1.0/Alpha)*rho[j]),(1.0/3.0));  }
   
   // smoothing mesh density function
  
    for (int j=1;j<N-1;j++) {
       rho_sm[j] = 0.25*( rh[j-1]+rh[j+1] )+0.5*rh[j]; }
       
   rho_sm[0] = 0.5*(rh[0]+rh[1]);
   
   rho_sm[N-1] = 0.5*(rh[N-1]+rh[N-2]);
  
 }   
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Mass = mass matrix /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void mass(double t, double *y, double **M)
   
   {
     int N=neq/2;

     vector<double> u(N);
     vector<double> x(N);

     for (int j=0;j<N;j++) {
		u[j] =y[j];
        x[j] =y[j+N]; } 

      /////// IC ///////

      double x0,u0,xNP1,uNP1;

      x0=0.0;
      u0=0.0;
      xNP1=1.0;
      uNP1=0.0;

      vector< vector<double> > M1(N,vector<double>(N));
      vector< vector<double> > M2(N,vector<double>(N));
      vector< vector<double> > M3(N,vector<double>(N));
      vector< vector<double> > M4(N,vector<double>(N));
      // M1
      for (int i=0;i<N;i++)  {
          M1[i][i]=1; }
      // M2
      M2[0][0]= -(u[1] - u0)/(x[1] - x0);

      M2[N-1][N-1]=- (uNP1 - u[N-2])/(xNP1 - x[N-2]);
         
      for (int i=0;i<N-1;i++) { 
		 M2[i][i]= - (u[i+1] - u[i-1])/(x[i+1] - x[i-1]); }

      // M3
      for (int i=0;i<N;i++) {
         for (int j=0;j<N;j++) {
           M3[i][j]=0; } }
       // M4
       M4[0][0]=-2;
       M4[0][1]=1;
       for (int i=1;i<N-1;i++) {
          M4[i][i]=-2;
          M4[i][i-1]=1;
          M4[i][i+1]=1; }
       M4[N-1][N-1]=-2;
       M4[N-1][N-2]=1;

        for (int i=0;i<N;i++) {
           for (int j=0;j<N;j++) {
               M[i][j]=M1[i][j];
               M[i][j+N]=M2[i][j];
               M[i+N][j]=M3[i][j];
               M[i+N][j+N]=M4[i][j]; } }
   }
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// gauss elimination //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void gauss(double *F, double **J, double *xn)

  { 
    vector< vector<double> > A(neq,vector<double>(neq+1)); // creating augmented matrix
 
    for (int j=0; j<neq; j++) {
       for (int k=0; k<neq; k++)  {
          A[j][k]=J[j][k]; } }

    for (int j=0; j<neq; j++)  {
       A[j][neq]=F[j]; }
       
    ///////// solving AX=B  ////////// 
    
    for (int i=0; i<neq; i++) {
        // Search for maximum in this column
       double maxEl = abs(A[i][i]);
       int maxRow = i;
       for (int k=i+1; k<neq; k++) {
          if (abs(A[k][i]) > maxEl) {
             maxEl = abs(A[k][i]);
             maxRow = k;  } }

    // Swap maximum row with current row (column by column)

       for (int k=i; k<neq+1;k++) {
       double tmp = A[maxRow][k];
       A[maxRow][k] = A[i][k];
       A[i][k] = tmp; }

    // Make all rows below this one 0 in current column

        for (int k=i+1; k<neq; k++) {
           double c = -A[k][i]/A[i][i];
           for (int j=i; j<neq+1; j++) {
              if (i==j) {
                 A[k][j] = 0; }
              else {
                 A[k][j] += c * A[i][j]; }
            } }
     }  // end i loop

    // Solve equation Ax=b for an upper triangular matrix A

    vector<double> x(neq);
   
    for (int i=neq-1; i>=0; i--) {
      x[i] = A[i][neq]/A[i][i];
      for (int k=i-1;k>=0; k--)  {
        A[k][neq] -= A[k][i] * x[i]; } }

    for (int k=0;k<neq; k++) {
      xn[k] = x[k]; }
	
  }

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////// rhs of the ode y'=L^{-1}g(t,y) /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void rhs(double t, double *y, double *f)

    {
      int N=neq/2;

      double ep, tau, x0,u0,xNP1,uNP1;

      ep=0.001;

      tau=0.1;

      vector<double> u(N);
      vector<double> x(N);

      for (int j=0;j<N;j++) {
        u[j] =y[j];
        x[j] =y[j+N]; } 
      /////// IC ///////
      x0=0.0;
      u0=0.0;
      xNP1=1.0;
      uNP1=0.0;
     ////////////////

      double *g = new double[neq];

      double dx;

      for (int i=1;i<(N-1);i++) {
         dx = x[i+1] - x[i-1];
         g[i] = (2.0*ep)/dx*( (u[i+1] - u[i])/(x[i+1] - x[i] ) - (u[i] - u[i-1])/(x[i] - x[i-1]) )- 0.5*(pow(u[i+1],2)-pow(u[i-1],2))/dx; }

         dx = x[1] - x0;    
    
         g[0] = (2.0*ep)/dx*((u[1] - u[0])/(x[1] - x[0]) - (u[0] - u0)/(x[0] - x0)) - 0.5*(pow(u[1],2) - pow(u0,2))/dx;

         dx = xNP1 - x[N-2];   
    
         g[N-1]=(2.0*ep)/dx*((uNP1 - u[N-1])/(xNP1 - x[N-1]) - (u[N-1] - u[N-2])/(x[N-1] - x[N-2]))/dx - 0.5*(pow(uNP1,2) - pow(u[N-2],2))/dx;
 
      double *rho_sm = new double[N]; 
   
      Rho(y,rho_sm); 

      double *v = new double[N];

      for (int i=1;i<(N-1);i++) {
         v[i] =( (rho_sm[i+1] + rho_sm[i])*(x[i+1] - x[i]) - (rho_sm[i] + rho_sm[i-1])*(x[i] - x[i-1]) ); }

         v[0] =  ( (rho_sm[1] + rho_sm[0])*(x[1] - x[0]) - (rho_sm[0] + rho_sm[0])*(x[0] - x0) );
    
         v[N-1] =( (rho_sm[N-1] + rho_sm[N-1])*(xNP1 - x[N-1]) - (rho_sm[N-1] + rho_sm[N-2])*(x[N-1] - x[N-2]) );

     for (int i=0;i<N;i++) {
        g[i+N]=-1.0/(2.0*tau)*v[i];}
    
      double **L = new double*[neq];

      for (int j=0;j<neq;j++) {
             L[j] = new double[neq]; }
      
      mass(t,y,L);  // calling mass matrix L(y) 

      double *w = new double[neq];

      gauss(g,L,w);  // linear solve

      for (int i=0;i<neq;i++) {
         f[i]=w[i]; }

      for (int i=0; i<neq; i++) {
       delete [] L[i]; }
  
      delete [] L;

      delete [] g;
    
      delete [] rho_sm;

      delete w;

      delete v;

   }

/////////////////////////////////////////////////////////////////////////////////////
////////////////////  function in the form F(x)=0  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    void fun(double t, double *p, double *u, double *Fn)

      {

        double* frh = new double[neq];

        rhs(t,p,frh);
   
        for (int j=0;j<neq;j++) {
           Fn[j]=p[j]-u[j]-dt*frh[j]; }

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
    
    for (int j=0;j<neq;j++) {
        xn[j]=xx[j]; }
     
    double *fx = new double[neq];

    fun(t,x,xold,fx);

    double *fxx = new double[neq];

    fun(t,xn,xold,fxx);

     for (int k=0;k<neq;k++) {
        J[k][i]=(fxx[k]-fx[k])/dx; }

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
          double accum = 0.0;
          for (int i = 0; i < n; ++i) {
            accum += w[i] * w[i]; }
          norm[0]=sqrt(accum);
        }
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Newton's Solver ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void step(double t, double *u, double *unew)

  {  
	  int max_iter=100;

      double tol=0.00000001;

	  int it_count=0;
	  
	  vector<double> p(neq);
      
      vector<double> z(neq);
        
      for (int j=0;j<neq;j++) {
        p[j]=u[j];  }
              
      double *w = new double[neq];   
        
      for (int j=0;j<neq;j++) {
        w[j]=u[j]; }
               
     // Newton's iteration starts here!

       for (int l=0;l<max_iter;l++)
            {
		    it_count=it_count+1;		
			     
             double *v = new double[neq];

             for (int j=0;j<neq;j++) {
                 v[j]=p[j]; }
             
             double *Fn = new double[neq];

             fun(t,v,w,Fn);   // calling function F(X)=0
             
             double **Jn = new double*[neq];

               for (int j=0;j<neq;j++) {
                  Jn[j] = new double[neq]; }
             
             jac(t,v,w,Jn);   // calling jacobian 

             double *x=new double[neq];

             gauss(Fn,Jn,x);  // calling linear solver for x=Jn/Fn
         
             for (int j=0;j<neq;j++) {
              z[j]=v[j]-x[j]; }
             
             double *dv=new double[neq];

             for (int j=0;j<neq;j++) {
               dv[j]=abs(z[j]-v[j]);}

             double *norm=new double[1]; 
              
             l2_norm(dv,neq,norm);
       
             if (norm[0]<tol) {
		     break; }
        
             p=z; 

             delete [] x;

             delete [] Fn;

             delete [] v;
          
             delete [] dv;
          
             delete [] norm;

             for (int i=0; i<neq; i++) {
             delete [] Jn[i]; }
             delete [] Jn;
       } // end newton loop l
       
            
     for (int j=0;j<neq;j++) {
      unew[j]=z[j]; } 
           
     delete [] w;
   
  } // end step

 
};


int main(int argc, char *argv[]) 

 {
  int start_s=clock(); 
  int order, nt,N;
  double *sol;

  if (argc != 4) {
    printf("usage: <executable> <order> <nt>  >  output_file\n");
    fflush(stdout);
    exit(1); } 
  else  {
    order = atoi(argv[1]); // order of method
    nt = atoi(argv[2]); // number of time steps
    N = atoi(argv[3]); }  // number of spatial mesh points
  
  int neq = 2*N;
  int ti = 0;
  int tf = 1;
  
  double dt = (double)(tf - ti)/double(nt); // compute dt
  
  // initialize ODE variable

  ImplicitOde *ode = new ImplicitOde(neq,nt,ti,tf,dt);

   double h = 1.0/(N+1);

   double *xint = new double[N];

   double *uintv = new double[N];

   for (int i=0;i<N;i++)
       {
       xint[i]=h*(i+1);
       uintv[i]=sin(2.0*M_PI*xint[i]) + 0.5*sin(M_PI*xint[i]);
       }

   sol = new double[neq];

   for (int i=0;i<N;i++)
       {
       sol[i]=uintv[i];
       sol[i+N]=xint[i];
       }
 
  // call ridc 

  ridc_be(ode, order, sol);
  
  for (int i = 0; i < neq; i++) {
   printf("%17.16f\n", sol[i]); }
   
  delete [] sol;
  delete [] xint;
  delete [] uintv;    

printf("Time taken: %.12fs\n", (double)(clock() - start_s)/CLOCKS_PER_SEC); 

 
  
}


