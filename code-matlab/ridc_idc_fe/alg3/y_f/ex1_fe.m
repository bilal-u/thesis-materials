
t_int=0;

t_final=1;

tspan=[t_int t_final];

% step size

dt=0.01;

% initial condition

y0=1;

% dimension of vector

addpath(genpath('/users/labnet5/gr8/buddin/research/project_MMPDE_RIDC/thesis/ridc/code/new_file/ridc_idc_fe/alg3/source_code1'));

N=length(y0);
  
format long

myf = @(t,x) f1(t,x);

u_st =fe(myf,y0,tspan,dt);

sol=u_st(:,end)
  
  