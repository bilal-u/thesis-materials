%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t_int = left end point of the time interval
% t_final = right end point of the time interval 
% p is the desired order of the RIDC methods 

clear all; close all; clc; 

t_int=0;

t_final=1;

tspan=[t_int t_final];

a=0.01;

p=2;

K=20;

N=19;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


           % Initial condition of the IVP %
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 1/(N+1); % excluding boundary values

xint = h*(1:N)';

uint = sin(2*pi*xint) + (1/2)*sin(pi*xint);

u0=[uint; xint];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % adding source codes to the current path %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%addpath(genpath('/home/buddin/resreach/mun/codes/ridc/new_file/ridc_idc_be/alg3/burgers_eq/source_code1'));

format long
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fuction f is defined by a separate scripts named f1.m

%myf = @(t,x) f  (t,x);

dt_st=[a;a/2;a/4;a/8];

n=length(dt_st);

u_st=zeros(2*N,n);

for m=1:n
    
dt=dt_st(m);    

u_st(:,m)=ridc(p,u0,tspan,dt,K);

end

E1=norm( abs(u_st(:,1)-u_st(:,2)) );

E2=norm( abs(u_st(:,2)-u_st(:,3)) );

E3=norm( abs(u_st(:,3)-u_st(:,4)) );


p1=log2(E1/E2);

p2=log2(E2/E3);

order=[p1;p2]
