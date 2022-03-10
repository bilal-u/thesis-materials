
tic;

clear all;

t_int=0;

t_final=1;

tspan=[t_int t_final];

% step size

dt=0.01;

p=1;

K=100;

% initial condition

y0=[0;1];

% dimension of vector

N=length(y0);

% adding source code to the current path

addpath(genpath('/home/buddin/resreach/thesis/ridc/code/new_file/ridc_idc_fe/alg3/My_f/source_code1'));

format long

myf = @(t,x) f1(t,x);

u_st =ridc(myf,p,y0,tspan,dt,K);

u_st(:,end)

toc;
  
   

