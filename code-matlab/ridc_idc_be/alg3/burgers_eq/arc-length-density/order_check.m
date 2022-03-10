
tic;

format long

clear all;

t_int=0;

t_final=1;

tspan=[t_int t_final];

% step size

a=0.001;

dt_ar=[a;a/2;a/4;a/8;a/16];

n=length(dt_ar);


p=2;

K=t_final/a;

N = 19;

h = 1/(N+1);

xinit = h*(1:N)';

uinit = sin(2*pi*xinit) + 0.5*sin(pi*xinit);

y0 = [uinit; xinit];

u_s=zeros(2*N,n);

for k=1:n
    
dt=dt_ar(k);

%%%%%%%%%%%%%%%%%%% RIDC %%%%%%%%%%%%%%%%%%%%%%

u_st =ridc(p,y0,tspan,dt,K);

u_s(:,k)=u_st(:,end);

end

E1=norm( abs(u_s(:,1)-u_s(:,2)) ,inf);

E2=norm( abs(u_s(:,2)-u_s(:,3)),inf );

E3=norm( abs(u_s(:,3)-u_s(:,4)),inf );

E4=norm( abs(u_s(:,4)-u_s(:,5)),inf );
 
% E5=norm( abs(u_s(:,5)-u_s(:,6)),inf );
% 
% E6=norm( abs(u_s(:,6)-u_s(:,7)),inf );

p1=log2(E1/E2);

p2=log2(E2/E3);

p3=log2(E3/E4);
 
% p4=log2(E4/E5);

% p5=log2(E5/E6);

orders=[0;p1;p2;p3]

errors=[E1;E2;E3;E4]

toc;
  
   

