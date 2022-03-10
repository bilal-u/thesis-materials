
function rh=ridc(f,p,y0,tspan,dt,K)

p1=p;

% K=p;

if p1==1
    
    rh=ridc1(f,y0,tspan,dt,K);
    
else

s=matrix(p1);  

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
    
   u(:,1,1)=y(:,j); 
    
    % prediction loop
    
    t(1)=t0;
    
%     for m=1:K
%         
%         t(m+1)=((j-1)*K+(m))*dt;
%           
%         u(:,m+1,1)= u(:,m,1) + dt*f(t(m),u(:,m,1));
%         
%     end


    
    % correction loop
    
    for l=1:M
        
       u(:,1,l+1)=u(:,1,l);
        
        for m=1:l
            
            t_0=((j-1)*K)*dt;
            
            s1=0; 
         
               for i=1:l+1
            
                  s1=s1+s(m,i,l)*f(t_0+(i-1)*dt,u(:,i,l));   
              
               end
            
         u(:,m+1,l+1)=u(:,m,l+1)+dt*( f(t(m),u(:,m,l+1))-f( t(m),u(:,m,l) ) ) +dt*s1; 
         
        end
        
        for m=l+1:K
          
             s2=0;
         
              for i=1:l+1
            
                 s2=s2+s(l,i,l)*f(t(m-l+i),u(:,m-l+i,l));   
         
              end
                          
            u(:,m+1,l+1)=u(:,m,l+1)+dt*( f(t(m),u(:,m,l+1))-f( t(m),u(:,m,l) ) ) + dt*s2;
            
        end    
        
    end
    
    y(:,j+1)=u(:,K+1,M+1);
    
    t0=t(K+1,1);
        
end 

rh=y;

end

end



   
    
    
