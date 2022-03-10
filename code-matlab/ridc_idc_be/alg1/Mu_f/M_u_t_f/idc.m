
function [orders,errors]=idc(a,p,s,T)

format long

M=p-1;

dt_ar=[a;a/2;a/4;a/8;a/16;a/32];

n=length(dt_ar);

% error storing

er=zeros(n,1)';

for k=1:n
    
dt=dt_ar(k);

W=T/dt; 

J=int16(W/M);

u0=[0;1];

u_st=zeros(2,J+1);

J=double(J);

Y=zeros(2,M+1,M+1);

u_st(:,1)=u0;

max_iter=100;

tol=1e-8;

for j=1:J 
    
    y0=u0;
    
   % prediction loop
    
   for m = 1:M
       
       Y(:,1,1)=y0;
       
       v=Y(:,m,1);
       
       t=((j-1)*M+m)*dt;
       
    % newton loop 
    
    for p=1:max_iter 
        
        Fn=fun_F(v,Y(:,m,1),dt,t);
        
        Jn=NumJacob(@fun_F,v,Y(:,m,1),dt,t); 
        
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
        
     Y(:,1,l+1)=Y(:,1,l); % initial guess
     
        for m=1:M
     
           v1=Y(:,m,l+1);
     
           i=m;
        
           q=l;
           
           t=((j-1)*M+m)*dt;
           
           t0=((j-1)*M+M)*dt;
           
           s0=inv(mass(t,Y(:,m+1,l)));
           
           s1=0;
           
           for n1=1:M+1
           
           s1=s1+s(m,n1)*inv(mass((t0-(M+1-n1)*dt),Y(:,n1,l)))*f((t0-(M+1-n1)*dt),Y(:,n1,l));
          
           end
% newton loop 
  
           for p=1:max_iter 
        
               Fn1=fun_corr(v1,Y(:,:,:),dt,i,q,s0,s1,t);
        
               Jn1=NumJacob(@fun_corr,v1,Y(:,:,:),dt,i,q,s0,s1,t);
                            
               z1=v1-Jn1\Fn1; 
        
                  if norm(abs(z1-v1),inf)<tol
               
                %   disp('converged')
                
                   break
                   
                  end
                  
             v1=z1;
             
           end
                 
       Y(:,m+1,l+1) = z1;
   
        end
        
   end
   
 u_st(:,j+1)= Y(:,M+1,M+1);  
   
 u0=u_st(:,j+1);
   
end 

u_ex=[sin(T);cos(T)];

er(k)=norm( abs(u_st(:,J+1)-u_ex),inf);

end

p1=log2(er(n-5)/er(n-4));

p2=log2(er(n-4)/er(n-3));

p3=log2(er(n-3)/er(n-2));

p4=log2(er(n-2)/er(n-1));

p5=log2(er(n-1)/er(n));

orders=[p1;p2;p3;p4;p5];
   
errors=er'; 

end