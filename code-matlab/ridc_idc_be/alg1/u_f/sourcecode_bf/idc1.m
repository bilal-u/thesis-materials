
function rh=idc1(f,u0,tspan,dt)

M=1;

N=length(u0);

T=tspan(1,2);

W=int16(T/dt); 

J=int16(W/M); 

y=zeros(N,J+1);

u=zeros(N,M+1,1);

y(:,1)=u0;

t=zeros(M+1,1);

J=double(J);

t0=tspan(1);

max_iter=100;

tol=1e-8;

for j=1:J 
    
    y0=u0;
    
    t(1)=t0;
       
    % prediction loop
    
    for m=1:M
        
        t(m+1)=((j-1)*M+m)*dt;
        
        u(:,1,1)=y0;
        
        v=u(:,m,1);
 
% start newton loop 

         for p=1:max_iter 
        
             Fn=fun_pred(v,u(:,m,1),t(m+1),dt,f);
        
             Jn=NumJacob(@fun_pred,v,u(:,m,1),t(m+1),dt,f);
        
             z=v-Jn\Fn; 
        
                  if norm(abs(z-v),inf)<tol
               
                 %  disp('converged')
                
                   break
                   
                  end
             v=z;
             
         end
         
 % end newton loop  
 
    u(:,m+1,1) = z;
                     
    end
    
   y(:,j+1)= u(:,M+1,1);  
   
   u0=y(:,j+1);
   
   t0=t(M+1,1);
    
end 

rh=y;

end

