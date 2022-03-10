
% objective of the program
              
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     this code written in matlab will find the numerical 
%     solution of a given IVP using RIDC methods contructed with
%     backward Euler integrators of orders from 1 to p in adaptive mesh.
%     Errors and the orders of convergenge of
%     methods are tabulated by a matrix. This code also create 
%     a command to create a table in latex.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
               % Physical problems %
                    
                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Consider the following nonlnear Burger's equation:

%    u_t=epsilon*u_xx - (u^2/2)_x,  x~[0,1], t~[0, 1]

%    initial condition: u(x,0)=sin(2*pi*x)+(1/2)*sin(pi*x);

%    Boundary condition: u(0,t)=u(1,t)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Input parameters %
                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t_int = left end point of the time interval
% t_final = right end point of the time interval 
% p is the desired order of the RIDC methods 

t_int=0;

t_final=.12;

tspan=[t_int t_final];

a=0.01;

p=3;

K=p;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


           % Initial condition of the IVP %
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y0 = initial condition
% N = the dimension of the vector y0 

N=19;

h = 1/(N+1);

xint = h*(1:N)';

uint = sin(2*pi*xint) + (1/2)*sin(pi*xint);

u0=[uint; xint];

W1=length(u0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % adding source codes to the current path %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('/home/buddin/resreach/mun/codes/ridc/new_file/ridc_idc_be/alg2/adaptive_mesh/ridc_be/source_code'));

format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       % function to the right side of the given ODE % 
         
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fuction f is defined by a separate scripts named f1.m

myf = @(t,x) f1  (t,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


          % exact solution from ode15s %


opts = odeset('RelTol',1e-18,'AbsTol',1e-16,'Mass',@mass,'MaxOrder',5);

sol= ode15s(@f1,tspan,u0,opts);

y = deval(sol,t_final);

% time step is halved % 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%  time step is halved and stored in an array in order 
%  to compute the order of the methods 

dt_ar=[a;a/2;a/4;a/8];

n=length(dt_ar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


er=zeros(n,1);


%while (m<=p)

% u_s=zeros(W1,n);
%  
% u=zeros(W1/2+2,n);
% x=zeros(W1/2+2,n);
% 
% Y=zeros(W1+4,1);

for k=1:n
    
   dt=dt_ar(k);   

   W=t_final/dt;

   u_st=ridc(myf,p,u0,tspan,dt,K);
   
   er(k)=norm(abs(y-u_st(:,end)), inf);
   
end

p1=log2(er(n-3)/er(n-2));

p2=log2(er(n-2)/er(n-1));

p3=log2(er(n-1)/er(n));
   
P=[p1;p2;p3]
% 
%    orders(2:n)=P(1:3);
%    
    errors=er;
   
%    orders_st(:,m)=orders(1:n);
%    
%    errors_st(:,m)=errors(1:n);
%    
%    m=m+1;
   
%    P=zeros(n-1,1);
%    
%    errors=zeros(n,1);
   
% craeting talex table 

% Er=errors_st;
% 
% Or=orders_st;
% 
% digits(4)
% 
% format short g
% 
% digits(4)
% 
% A=zeros(n,2*p);
% 
% A(:,1)=Er(:,1);
% A(:,2)=Or(:,1);
% 
% A(:,3)=Er(:,2);
% A(:,4)=Or(:,2);
% 
% A(:,5)=Er(:,3);
% A(:,6)=Or(:,3);
% 
% A(:,7)=Er(:,4);
% A(:,8)=Or(:,4);
% 
% Vector= A
% 
% digits(4)
%  
% result=latex(sym(vpa(A)))
%  
%  
