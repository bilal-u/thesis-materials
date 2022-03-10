
% objective of the program
              
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     this code written in matlab will find the numerical 
%     solution of a given IVP using IDC methods contructed with
%     forward Euler integrators of orders from 1 to p.
%     Errors and the orders of convergenge of
%     methods are tabulated by a matrix. This code also create 
%     a command to create a table in latex.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
               % Physical problems %
                    
                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Consider the following linear advection-diffusion equation:

%    u_t=u_x + epsilon*u_xx,  x~[0,1], t~[0, 1.2]

%    initial condition: u(x,0)=2 + sin(2*pi*x);

%    with periodic Boundary conditions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


             % Input parameters %
                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t_int = left end point of the time interval
% t_final = right end point of the time interval 
% p is the desired order of the IDC methods 
% a is the step size 

t_int=0;

t_final=1.2;

tspan=[t_int t_final];

a=0.01;

p=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


           % Initial condition of the IVP %
                 
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y0 = initial condition
% N = the dimension of the vector y0 

N=10;

x=linspace(0,1,N)';

y0=2+sin(2*pi*x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       % function to the right side of the given ODE % 
         
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fuction f is defined by a separate scripts named f1.m

myf = @(t,x) f6  (t,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

        % adding source codes to the current path %
        
        
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/users/labnet5/gr8/buddin/research/project_MMPDE_RIDC/thesis/ridc/code/new_file/ridc_idc_be/uniform_mesh/idc_be/u_f/sourcecode_bf'));

format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      % time step is halved % 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%  time step is halved and stored in an array in order 
%  to compute the order of the methods 

dt_ar=[a;a/2;a/4;a/8;a/16];

n=length(dt_ar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      % Initializztion of vectors %
      
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_st=zeros(N,n);

orders=zeros(n-1,1);

errors=zeros(n-1,1);

orders_st=zeros(n-1,p);

errors_st=zeros(n-1,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


m=1; % counter

while (m<=p) % loop runs from 1 to p 
    

   for k=1:n
       
   % step size     
    
   dt=dt_ar(k);   
   
   % calling idc

   u_s=idc(myf,m,y0,tspan,dt);
   
   u_st(:,k)=u_s(:,end);
   
   end
   
   % infinity norm of the relative errors
   
   
   E1=norm( abs(u_st(:,1)-u_st(:,2)) ,inf);

   E2=norm( abs(u_st(:,2)-u_st(:,3)),inf );

   E3=norm( abs(u_st(:,3)-u_st(:,4)),inf );

   E4=norm( abs(u_st(:,4)-u_st(:,5)),inf );

   p1=log2(E1/E2);

   p2=log2(E2/E3);

   p3=log2(E3/E4);

   % orders and errors storing

   orders=[0;p1;p2;p3];

   errors=[E1;E2;E3;E4];

   orders_st(:,m)=orders(1:n-1);
   
   errors_st(:,m)=errors(1:n-1);
   
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

%digits(4)

format short g

%digits(4)

A=zeros(n-1,2*p);

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

%result=latex(sym(vpa(A)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %  The end % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
