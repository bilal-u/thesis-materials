
function rh=fe(f,y0,tspan,dt)

N=length(y0);

T=tspan(1,2);

t1=tspan(1):dt:T;

l=length(t1);

u=zeros(N,l);

u(:,1)=y0;

% t0=0;

   
    for m=1:l-1
        
       % t(m+1)=t(m)+m*dt;
        
        u(:,m+1)= u(:,m) + dt*f(t1(m),u(:,m));
    end
   
rh=u;

end