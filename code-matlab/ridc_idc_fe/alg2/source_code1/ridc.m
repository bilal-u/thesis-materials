
function rh=ridc(f,p,y0,tspan,dt)

i=p;

K=p;

if i==1
    
    rh=ridc1(f,y0,tspan,dt,K);
    
else

s=integration_matrix(i);  

N=length(y0);

M=p-1;

T=tspan(1,2);

W=int16(T/dt); 

J=int16(W/K); 

y=zeros(N,J+1);

u=zeros(N,K+1,p);

y(:,1)=y0;

t=zeros(K+1,1);

J=double(J);

M=double(M);

t0=tspan(1);

for j=1:J 
    
   u(:,1,1)=y(:,j); % need to be generelized
    
    % prediction loop
    
    t(1)=t0;
    
    for m=1:K
        
        t(m+1)=((j-1)*K+(m))*dt;
          
        u(:,m+1,1)= u(:,m,1) + dt*f(t(m),u(:,m,1));
        
    end
    
    % correction loop
    
    for l=1:M
        
       u(:,1,l+1)=u(:,1,l);
        
        for m=1:M
            
            t_0=((j-1)*K)*dt;
            
            s1=0; 
         
               for n1=1:M+1
            
                  s1=s1+s(m,n1)*f(t_0+(n1-1)*dt,u(:,n1,l));   
              
               end
            
         u(:,m+1,l+1)=u(:,m,l+1)+dt*( f(t(m),u(:,m,l+1))-f( t(m),u(:,m,l) ) ) + M*dt*s1; 
         
        end
        
        for m=M+1:K
          
             s2=0;
         
              for n1=1:M+1
            
                 s2=s2+s(M,n1)*f(t(m-M+n1),u(:,m-M+n1,l));   
         
              end
                          
            u(:,m+1,l+1)=u(:,m,l+1)+dt*( f(t(m),u(:,m,l+1))-f( t(m),u(:,m,l) ) ) + M*dt*s2;
            
        end    
        
    end
    
    y(:,j+1)=u(:,K+1,M+1);
    
    t0=t(K+1,1);
    
end 

rh=y;

end

end



   
    
    