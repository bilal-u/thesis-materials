

function rh=ridc(f,p,y0,tspan,dt,K)

i=p;

%K=i;

  if i==1
    
    rh=ridc1(f,y0,tspan,dt,K);
    
  else

s=matrix(i);  

% initial condition & guess

N=length(y0);

M=p-1;

T=tspan(1,2);

W=int16(T/dt); 

J=int16(W/K); 

t=zeros(K+1,1);

t0=tspan(1);

u_st=zeros(N,J+1);

Y=zeros(N,K+1,M+1);

u_st(:,1)=y0;

max_iter=100;

tol=1e-8;

for j=1:J 
    
   t(1)=t0;
    
   % prediction loop
    
   for m = 1:K
       
       Y(:,1,1)=y0;
       
       v=Y(:,m,1);
       
       t(m+1)=((j-1)*K+m)*dt;
       
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
   
     % correction loop
    
    for l=1:M
        
     Y(:,1,l+1)=Y(:,1,l); 
     
        for m=1:l
     
           v1=Y(:,m,l+1);
     
           i=m;
        
           q=l;
           
           s0=inv(mass(t(m+1),Y(:,m+1,l))) *f(t(m+1),Y(:,m+1,l));
           
           s1=0;
           
           t_0=((j-1)*K)*dt;
           
           for n1=1:l+1
           
           s1=s1+s(m,i,l)*inv(mass(t_0+(n1-1)*dt,Y(:,n1,l)))*f(t_0+(n1-1)*dt,Y(:,n1,l));
          
           end
          
% newton loop 
  
           for p=1:max_iter 
        
               Fn1=fun_corr(v1,Y(:,:,:),dt,t(m+1),i,q,s0,s1,f);
        
               Jn1=NumJacob(@fun_corr,v1,Y(:,:,:),dt,t(m+1),i,q,s0,s1,f);
                            
               z1=v1-Jn1\Fn1; 
        
                  if norm(abs(z1-v1),inf)<tol
               
                 %  disp('converged')
                
                   break
                   
                  end
                  
             v1=z1;
             
           end
           
       Y(:,m+1,l+1) = z1;
   
        end
        
       for m=l+1:K
     
           v2=Y(:,m,l+1);
     
           i=m;
        
           q=l;
           
           s2=inv(mass(t(m+1),Y(:,m+1,l))) * f(t(m+1),Y(:,m+1,l));
           
           s3=0;
           
           for n1=1:l+1
           
           s3=s3+s(l,i,l)*inv(mass(t(m-l+n1),Y(:,m-l+n1,l)))*f(t(m-l+n1),Y(:,m-l+n1,l));
           
           end
         
% newton loop 
  
           for p=1:max_iter 
               
               Fn2=fun_corr(v2,Y(:,:,:),dt,t(m+1),i,q,s2,s3,f);
        
               Jn2=NumJacob(@fun_corr,v2,Y(:,:,:),dt,t(m+1),i,q,s2,s3,f);
                            
               z2=v2-Jn2\Fn2; 
        
                  if norm(abs(z2-v2),inf)<tol
               
                  % disp('converged')
                
                   break
                   
                  end
                  
             v2=z2;
             
           end
           
       Y(:,m+1,l+1) = z2;
   
       end   
        
   end
   
 u_st(:,j+1)= Y(:,K+1,M+1);  
   
 y0=Y(:,K+1,M+1);
 
 t0=t(K+1,1);
   
end 

rh=u_st;

 end

end





