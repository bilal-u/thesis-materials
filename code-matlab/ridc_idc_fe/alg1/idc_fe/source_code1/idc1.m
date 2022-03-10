
function rh=idc1(f,y0,tspan,dt)

M=1;

N=length(y0);

T=tspan(1,2);

W=int16(T/dt); 

J=int16(W/M); 

y=zeros(N,J+1);

u=zeros(N,M+1,1);

y(:,1)=y0;

t0=tspan(1);

t=zeros(M+1,1);

J=double(J);

for j=1:J 
    
    u(:,1,1)=y(:,j); 
    
    % prediction loop
    
    t(1)=t0;
    
    for m=1:M
        
        t(m+1)=((j-1)*M+(m-1))*dt;
        
        u(:,m+1,1)= u(:,m,1) + dt*f(t(m),u(:,m,1));
    end
    
   
    y(:,j+1)=u(:,M+1,1);
    
   
    t0=t(M+1,1);
    
end
   
rh=y;

end

