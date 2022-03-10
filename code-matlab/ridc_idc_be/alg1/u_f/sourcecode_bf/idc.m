
function rh=idc(f,p,u0,tspan,dt)

i=p;

%K=i;

if i==1
    
    rh=idc1(f,u0,tspan,dt);
    
else

s=integration_matrix(i); 

N=length(u0);

M=p-1;

T=tspan(1,2);

W=int16(T/dt); 

J=int16(W/M); 

y=zeros(N,J+1);

u=zeros(N,M+1,p);

y(:,1)=u0;

t=zeros(M+1,1);

J=double(J);

M=double(M);

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
    
     % correction loop
    
    for l=1:M
        
     u(:,1,l+1)=u(:,1,l); % initial guess
     
        for i=1:M
            
        t_0=((j-1)*M)*dt;      
     
        v=u(:,i,l+1);
     
        m=i;
        
        q=l; 
        
         s1=0;
        
          for n1=1:M+1
           
           s1=s1+s(m,n1)*f(t_0+(n1-1)*dt,u(:,n1,l));
          
           end
 
% newton loop 
  
           for p=1:max_iter 
               
               Fn=fun_corr(v,u(:,:,:),t(m+1),dt,s1,m,q,f,M);
        
               Jn=NumJacob(@fun_corr,v,u(:,:,:),t(m+1),dt,s1,m,q,f,M);
        
               z=v-Jn\Fn; 
        
                  if norm(abs(z-v),inf)<tol
               
                %   disp('converged')
                
                   break
                   
                  end
             v=z;
             
           end
           
       u(:,i+1,l+1) = z;
   
        end
      
    end
   
   y(:,j+1)= u(:,M+1,M+1);  
   
   u0=u(:,M+1,M+1);
   
   t0=t(M+1,1);
    
end 

rh=y;

end

end

