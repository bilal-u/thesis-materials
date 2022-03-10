

%clear; clc;

format long

T=1;

%dt=0.01;

J=101;

dt=0.01;

Mass=[4 -1; 
     -1 4];

y=zeros(2,J);

h=zeros(2,1);

u0=[0;1]; % IC  

y(:,1)=u0;

f(y(:,1))


for j=2:J
    
    h=f(y(:,j-1))
    
    y(:,j)=y(:,j-1) + dt*(Mass\h);%f(y(:,j-1)));
    
end 

y(:,J)

%e1=norm(y(:,J),inf);

%e2=norm(y(:,J),inf);