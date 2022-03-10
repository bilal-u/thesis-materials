
function rh=ridc1(y0,tspan,dt,K)

N=length(y0);

T=tspan(1,2);

W=int16(T/dt); 

J=int16(W/K); 

y=zeros(N,J+1);

u=zeros(N,K+1,1);

y(:,1)=y0;

t0=tspan(1);

t=zeros(K+1,1);

J=double(J);

max_iter=100;

tol=1e-8;

it_count=zeros(K,1);

for j=1:J 
    
    u(:,1,1)=y(:,j);  % copying initial solution to each group
    
    % prediction loop
    
    t(1)=t0;  % time initialization  in each group 
    
   
  for m=1:K
        
        v=u(:,m,1);
        
        t(m+1)=t(1)+m*dt;
        
        it=0;

      for p=1:max_iter 
        
        Fn=fun_F(v,u(:,m,1),dt);
        
        Jn=jacobFD(@fun_F,v,u(:,m,1),dt); 
       
        z=v-Jn\Fn; 
        
        it=it+1;
       
           if norm(abs(z-v))<tol
          
           %  disp('converged')
             break
             
           end
           
        v=z;
        
      end
    
    u(:,m+1,1) = z;
    
    it_count(m)=it;
    
  end    
  
 y(:,j+1)=u(:,K+1,1); 

 t0=t(K+1,1);
    
end
   
rh=y(:,end);

%it_count

end

