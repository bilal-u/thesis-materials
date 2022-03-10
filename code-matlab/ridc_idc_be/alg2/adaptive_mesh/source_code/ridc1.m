

function rh=ridc1(f,y0,tspan,dt,K)

N=length(y0);

T=tspan(1,2);

W=int16(T/dt); 

J=int16(W/K); 

t=zeros(K+1,1);

t0=tspan(1);

u_st=zeros(N,J+1);

Y=zeros(N,K+1,1);

u_st(:,1)=y0;

max_iter=100;

tol=1e-8;

for j=1:J 
    
   t(1)=t0;
    
   % prediction loop
    
   for m = 1:K
       
       Y(:,1,1)=y0;
       
       v=Y(:,m,1);
       
       t(m+1)=((j-1)*K+(m))*dt;
       
    % newton loop 
    
    for p=1:max_iter 
        
        Fn=fun_F(v,Y(:,m,1),dt,t(m+1),f);
        
        Jn=NumJacob(@fun_F,v,Y(:,m,1),dt,t(m+1),f); 
        
        z=v-Jn\Fn; 
        
         if norm(abs(z-v),inf)<tol
             
            % disp('converged')
             
             break
             
         end
         
        v=z;
        
    end
    
    Y(:,m+1,1) = z;
    
  end    
   
  
 u_st(:,j+1)= Y(:,K+1,1);  
   
 y0=u_st(:,j+1);
 
 t0=t(K+1,1);
   
end 

rh=u_st;

end

