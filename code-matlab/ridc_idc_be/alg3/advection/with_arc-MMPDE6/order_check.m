tic;

format long

clear all;

t_int=0.0;

t_final=1.0;

tspan=[t_int t_final];

%dt=0.01;

p=1;

K=20;

% step size

a=0.01;

dt_ar=[a;a/2;a/4;a/8;a/16];

n=length(dt_ar);


%K=t_final/a;

N = 39;

L = 10;

x0 = L/2;

sg = 1/2;

h = L/(N+1);

xint = h*(1:N)';

uint = exp(-((xint-x0).^2)/sg); %sin(2*pi*xint) + (1/2)*sin(pi*xint);

%u0=[uint; xint];

%uinit = sin(2*pi*xinit) + 0.5*sin(pi*xinit);

y0 = [uint; xint];

u_s=zeros(2*N,n);

for k=1:n
    
dt=dt_ar(k);

%%%%%%%%%%%%%%%%%%% RIDC %%%%%%%%%%%%%%%%%%%%%%

u_st =ridc(p,y0,tspan,dt,K);

u_s(:,k)=u_st(:,end);

end

E1=norm( abs(u_s(:,1)-u_s(:,2)));

E2=norm( abs(u_s(:,2)-u_s(:,3)));

E3=norm( abs(u_s(:,3)-u_s(:,4)));

E4=norm( abs(u_s(:,4)-u_s(:,5)));
 
% E5=norm( abs(u_s(:,5)-u_s(:,6)),inf );
% 
% E6=norm( abs(u_s(:,6)-u_s(:,7)),inf );

p1=log2(E1/E2);

p2=log2(E2/E3);

p3=log2(E3/E4);
 
% p4=log2(E4/E5);

% p5=log2(E5/E6);

order=[0;p1;p2;p3]

error=[E1;E2;E3;E4]

toc;
