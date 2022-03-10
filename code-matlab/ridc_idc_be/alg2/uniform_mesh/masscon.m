
% objective of the program
              
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     this code written in matlab will find the numerical 
%     solution of a given IVP using RIDC methods contructed with
%     backward Euler integrators of orders from 1 to p.
%     Errors and the orders of convergenge of
%     methods are tabulated by a matrix. This code also create 
%     a command to create a table in latex.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
               % Physical problems %
                    
                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Consider the IVP of the form :

%    Ly^{'}(t)=g(t,y), \; t\in [0,1.2]

%L=[4  -1 \\
%  -1  4 ]
%g(t,y)=(y_1+4y_2,-4y_1-y_2)^{t}.


%The exact solution is given by $y(t)=(\sin(t),\cos(t))^T$.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


             % Input parameters %
                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t_int = left end point of the time interval
% t_final = right end point of the time interval 
% p is the desired order of the RIDC methods 
% a is the step size 

t_int=0;

t_final=1.2;

tspan=[t_int t_final];

a=0.02;

p=4;

K=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


           % Initial condition of the IVP %
                 
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y0 = initial condition
% N = the dimension of the vector y0 

y0=[0;1];

N=length(y0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


         % function to the right side of the given ODE % 
         
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fuction f is defined by a separate scripts named f1.m

myf = @(t,x) massconfun  (t,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             
        % True solution if available %
                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The true solution of this IVP is :

%        y(t)=[cos(t);sin(t)];

u_ex=zeros(N,1);

u_ex(:,1)=[sin(t_final);cos(t_final)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


        % adding source codes to the current path %
        
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/users/labnet5/gr8/buddin/research/ridc/new_file/ridc_idc_be/alg2/uniform_mesh/u_f/sourcecode_bf'));

format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      % time step is halved % 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%  time step is halved and stored in an array in order 
%  to compute the order of the methods 

dt_ar=[a;a/2;a/4;a/8];

n=length(dt_ar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      % Initializztion of vectors %
      
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

er=zeros(n,1);

orders=zeros(n,1);

errors=zeros(n,1);

orders_st=zeros(n,p);

errors_st=zeros(n,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


m=1; % counter

while (m<=p) % loop runs from 1 to p 
    

   for k=1:n
       
   % step size     
    
   dt=dt_ar(k);   
   
   % calling ridc

   u_st=ridc(myf,m,y0,tspan,dt,K);
   
   % infinity norm error
   
   er(k)=norm( abs(u_st(:,end)-u_ex),inf);

   end
   
   % otders calculation
  
   p1=log2(er(n-3)/er(n-2));

   p2=log2(er(n-2)/er(n-1));

   p3=log2(er(n-1)/er(n));
   
   P=[p1;p2;p3];
   
   % orders and errors storing

   orders(2:n)=P(1:3);
   
   errors=er;
   
   orders_st(:,m)=orders(1:n);
   
   errors_st(:,m)=errors(1:n);
   
   %  moving to the (p+1)th order from p^th order method
   
   m=m+1;
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       % tabulating errors and the orders of the methods % 
       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    errors and orders of the corresponding p^th order methods
%    are merged in a single matrix.    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Er=errors_st;

Or=orders_st;

digits(4)

format short g

digits(4)

A=zeros(n,2*p);

A(:,1)=Er(:,1);
A(:,2)=Or(:,1);

A(:,3)=Er(:,2);
A(:,4)=Or(:,2);

A(:,5)=Er(:,3);
A(:,6)=Or(:,3);

A(:,7)=Er(:,4);
A(:,8)=Or(:,4);

Vector= A; 

disp(Vector)


digits(4)
 
result=latex(sym(vpa(A)))
 
 
