
function rh=ridc(f,p,u0,tspan,dt,K)

m=p;

if m==1
    
    rh=ridc1(f,u0,tspan,dt,K);
    
else

s=integration_matrix(m); 

N=length(u0);

M=p-1;

T=tspan(1,2);

W=int16(T/dt); 

J=int16(W/K); 

y=zeros(N,J+1);

u=zeros(N,K+1,M+1);

y(:,1)=u0;

t=zeros(K+1,1);

J=double(J);

M=double(M);

t0=tspan(1);

max_iter=100;

tol=1e-8;

for j=1:J 
    
    u(:,1,1)=y(:,j);
    
    t(1)=t0;
       
    % prediction loop
    
    for m=1:K
        
        t(m+1)=((j-1)*K+(m-1))*dt;
        
        v=u(:,m,1);
 
% start newton loop 

         for p=1:max_iter 
        
             Fn=fun_pred(v,u(:,m,1),t(m+1),dt,f);
        
             Jn=jacobFD(@fun_pred,v,u(:,m,1),t(m+1),dt,f);
        
             z=v-Jn\Fn; 
        
                  if norm(abs(z-v))<tol
               
                     break
                   
                  end
             v=z;
             
         end
         
 % end newton loop  
 
    u(:,m+1,1) = z;
                     
    end
    
     % correction loop
    
    for l=1:M
           
      u(:,1,l+1)=u(:,1,l); 
     
        for m=1:l
       
            v=u(:,m,l+1);
               
            s1=0;
            
           for n1=1:l+1
           
              s1=s1+s(m,n1,l)*f(t(1),u(:,n1,l));
                    
           end
           
           %sm=u(:,m,l+1)+l*dt*s1-dt*f(t(1),u(:,m+1,l));
          
           
% newton loop 
  
           for p=1:max_iter 
               
               Fn=fun_corr(v,t(m+1),dt,s1,f);
        
              Jn=jacobFD(@fun_corr,v,t(m+1),dt,s1,f);
        
                         
               z=v-Jn\Fn; 
        
                  if norm(abs(z-v))<tol
                 
                     break
                   
                  end
             v=z;
             
           end
           
       u(:,m+1,l+1) = z;
   
        end
   
      for m=l+1:K
     
         v=u(:,m,l+1);
     
           s2=0; 
           
           for n1=1:l+1
           
               s2=s2+s(l,n1,l)*f(t(m-l+n1),u(:,m-l+n1,l));
          
           end
           
           %sm2=u(:,m,l+1)+l*dt*s2-dt*f(t(1),u(:,m+1,l));
           
         
% newton loop 
  
           for p=1:max_iter 
        
               Fn=fun_corr(v,t(m+1),dt,s2,f);
        
               Jn=jacobFD(@fun_corr,v,t(m+1),dt,s2,f);
        
                   
               z=v-Jn\Fn; 
        
                  if norm(abs(z-v))<tol
               
                     break
                   
                  end
             v=z;
             
           end
           
       u(:,m+1,l+1) = z;
   
       end
    
    end
   
   y(:,j+1)= u(:,K+1,M+1);  
   
   t0=t(K+1,1);
    
end 

rh=y;

end

end

