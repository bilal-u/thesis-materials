%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t_int = left end point of the time interval
% t_final = right end point of the time interval 
% p is the desired order of the RIDC methods 

clear all; clc; 

t_int=0.0;

t_final=1.0;

tspan=[t_int t_final];

dt=0.01;

p=1;

K=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           % Initial condition of the IVP %
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y0 = initial condition
% N = the dimension of the vector y0 

N=19;

h = 1/(N+1);

%x=zeros(N,1);

x=h*(1:N)';

uint = sin(2*pi*x) + (1/2)*sin(pi*x);

u0=[uint; x];

 % adding source codes to the current path %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('/users/labnet5/gr8/buddin/code_regular_update/ridc/new_file/ridc_idc_be/alg3/burgers_eq/fixed/source_code1'));

format long 
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fuction f is defined by a separate scripts named f1.m

myf = @(t,x) f (t,x);

u_st=ridc(myf,p,u0,tspan,dt,K);

sol=u_st(:,end);

 
