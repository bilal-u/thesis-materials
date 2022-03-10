
function rh=ridc(p,y0,tspan,dt,K)

p1=p;

 if p1==1
    
    rh=ridc1(y0,tspan,dt,K);
    
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

max_iter=100;

tol=1e-8;

t0=tspan(1);

for j=1:J
    
    u(:,1,1)=y(:,j);  % copying initial solution to each group
    
    % prediction loop
    
    t(1)=t0;          % time initialization  in each group 

      
    for m=1:K
        
         v=u(:,m,1);
        
        t(m+1)=t(1)+m*dt;

        for p=1:max_iter 
        
        Fn=fun_F(v,u(:,m,1),dt);
        
        Jn=jacobFD(@fun_F,v,u(:,m,1),dt);
        
        z=v-Jn\Fn; 
        
           if norm(abs(z-v))<tol
             
          %   disp('converged')
             
            break
                   
           end
           
            if p==max_iter
                disp('not converged')
            end   
         
        v=z;
        
       end
    
    u(:,m+1,1) = z;
    
    end    

    
    % correction loop
    
    for l=1:M
        
       u(:,1,l+1)=u(:,1,l);

       for m=1:l
     
           v1=u(:,m,l+1);
     
           i=m;
        
           q=l;
           
           s0=inv(mass(u(:,m+1,l)));
           
           s1=0;
           
            for n1=1:l+1
           
           s1=s1+s(m,n1,l)*inv(mass(u(:,n1,l)))*f(u(:,n1,l));
          
           end


       % newton loop 
  
           for p=1:max_iter 
        
               Fn1=fun_corr(v1,u(:,:,:),dt,i,q,s0,s1);
        
               Jn1=jacobFD(@fun_corr,v1,u(:,:,:),dt,i,q,s0,s1);
                            
               z1=v1-Jn1\Fn1; 
        
                  if norm(abs(z1-v1))<tol
               
                 %  disp('converged')
                
                   break
                   
                  end
                  
             v1=z1;
             
           end
           
       u(:,m+1,l+1) = z1;
   
        end
       
   for m=l+1:K
     
           v2=u(:,m,l+1);
     
           i=m;
        
           q=l;
           
           s2=inv(mass(u(:,m+1,l)));
           
           s3=0;
           
           for n1=1:l+1
           
           s3=s3+s(l,n1,l)*inv(mass(u(:,m-l+n1,l)))*f(u(:,m-l+n1,l));
           
           end
          
% newton loop 
  
           for p=1:max_iter 
        
               Fn2=fun_corr(v2,u(:,:,:),dt,i,q,s2,s3);
        
               Jn2=jacobFD(@fun_corr,v2,u(:,:,:),dt,i,q,s2,s3);
                            
               z2=v2-Jn2\Fn2; 
        
                  if norm(abs(z2-v2))<tol
               
                  % disp('converged')
                
                   break
                   
                  end
                  
             v2=z2;
             
           end
           
       u(:,m+1,l+1) = z2;
   
    end      
        
   end
   
    
    y(:,j+1)=u(:,K+1,M+1); % updating solution for the next group
    
    t0=t(K+1,1);           % updating time for the next group
    
end 

rh=y(:,end);

end

end



   
    
    
