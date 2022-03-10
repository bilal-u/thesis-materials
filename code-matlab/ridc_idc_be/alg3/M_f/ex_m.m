
tic;

t_int=0;

t_final=1;

tspan=[t_int t_final];

% step size

a=0.01;

p=1;

K=20;

% time step is halved 

dt_ar=[a;a/2;a/4;a/8];

n=length(dt_ar);

% error storing

er=zeros(n,1);

% initial condition

y0=[0;1];

% dimension of vector

N=length(y0);

% highest expected order of RIDC

errors=zeros(n,1);

 orders_st=zeros(n,p);
 
 errors_st=zeros(n,p);

% adding source code to the current path

addpath(genpath('/home/buddin/resreach/thesis/ridc/code/new_file/ridc_idc_fe/alg3/My_f/source_code1'));

format long

% exact solution at the right most point

u_ex=zeros(N,1);

u_ex(:,1)=[sin(t_final);cos(t_final)];

myf = @(t,x) f1(t,x);

 m=1;

while (m<=p)
    

   for k=1:n
    
   dt=dt_ar(k);   

   u_st =ridc(myf,p,y0,tspan,dt,K);
  
   er(k)=norm( abs(u_st(:,end)-u_ex),inf);

   end

    p1=log2(er(n-3)/er(n-2));
 
    p2=log2(er(n-2)/er(n-1));

    p3=log2(er(n-1)/er(n))
 
 toc;
   
   orders=[0;p1;p2;p3]

   errors=er;
   
    orders_st(:,m)=orders(1:n);
    
    errors_st(:,m)=errors(1:n);
   
 
% craeting talex table 

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
 
 Vector= A
 
 digits(4)
  
 %result=latex(sym(vpa(A)))
