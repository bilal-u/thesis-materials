
function rh=ridc(f,p,y0,tspan,dt,K)

p1=p;

 if p1==1
    
    rh=ridc1(f,y0,tspan,dt,K);
    
 else

s=integration_matrix(p1);  

N=length(y0);

M=p-1;

T=tspan(1,2);

W=int16(T/dt); 

J=int16(W/K);       % number of group is 1

y=zeros(N,J+1);     % storing solution groupwise

u=zeros(N,K+1,M+1); % temporary solution

y(:,1)=y0;

t=zeros(K+1,1);

J=double(J);

t0=tspan(1);

for j=1:J
    
    u(:,1,1)=y(:,j);  % copying initial solution to each group
    
    % prediction loop
    
    t(1)=t0;          % time initialization  in each group 
    
    for m=1:K
        
        t(m+1)=t(1)+m*dt;
        
        u(:,m+1,1)= u(:,m,1) + dt*f(t(m),u(:,m,1));
        
    end
    
    % correction loop
    
    for l=1:M
        
       u(:,1,l+1)=u(:,1,l);
        
        for m=1:l
            
            s1=0; 
         
               for i=1:l+1
            
                  s1=s1+s(m,i,l)*f(t(1)+(i-1)*dt,u(:,i,l));   
              
               end
            
         u(:,m+1,l+1)=u(:,m,l+1)+dt*( f(t(m),u(:,m,l+1))-f( t(m),u(:,m,l) ) ) +...
             l*dt*s1; 
         
        end
        
        for m=l+1:K
          
             s2=0;
         
              for i=1:l+1
            
                 s2=s2+s(l,i,l)*f(t(m-l+i),u(:,m-l+i,l));   
         
              end
                          
            u(:,m+1,l+1)=u(:,m,l+1)+dt*( f(t(m),u(:,m,l+1))-f( t(m),u(:,m,l) ) ) +...
                l*dt*s2;
            
        end    
        
    end
    
    y(:,j+1)=u(:,K+1,M+1); % updating solution for the next group
    
    t0=t(K+1,1);           % updating time for the next group
    
end 

rh=y;

end

end



   
    
    
