
t_int=0;

t_final=1;

tspan=[t_int t_final];

% step size

dt=0.01;

p=4;

K=100;

% initial condition

y0=1;

% dimension of vector

N=length(y0);
  
% adding source code to the current path

addpath(genpath('/users/labnet5/gr8/buddin/research/project_MMPDE_RIDC/thesis/ridc/code/new_file/ridc_idc_fe/test_alg3/source_code1'));

format long

myf = @(t,x) f2(t,x);

u_st =ridc(myf,p,y0,tspan,dt,K);

sol=u_st(:,end)
  
  