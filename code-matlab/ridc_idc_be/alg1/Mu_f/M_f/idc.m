
function [orders,errors]=idc(a,p,s,T)

format long

M=p-1;

% time step is halved 

dt_ar=[a;a/2;a/4;a/8;a/16;a/32];

%dt_ar=[a;a/2;a/4;a/8;a/16;a/32;a/64;a/128];

n=length(dt_ar);

% error storing

er=zeros(n,1);

u_st=zeros(2,n);

% mass matrix

Mass=[4 -1; 
     -1 4];

for k=1:n

dt=dt_ar(k);   

N=T/dt; 

J=int16(N/M);

u_ex=zeros(2,M+1,M+1);

u=zeros(2,M+1,M+1);

y=zeros(2,J+1);

u0=[0;1]; % IC  

y(:,1)=u0;

max_iter=100;

tol=1e-8;

for j=1:J 
    
    y0=u0;
       
    % prediction loop
    
    for m=1:M
        
        u(:,1,1)=y0;
        
       % t=((j-1)*M+(m))*dt;
        
        v=u(:,m,1);
 
% start newton loop 

         for p=1:max_iter 
        
             Fn=fun_pred(v,u(:,m,1),dt,Mass);
        
             %Jn=jac(v,dt); % analytic jacobian
             
             Jn=NumJacob(@fun_pred,v,u(:,m,1),dt,Mass); % numerical jacobian
        
             z=v-Jn\Fn; 
        
                   if norm(abs(z-v),inf)<tol
%                
%                    disp('converged')
                 
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
     
        v=u(:,i,l+1);
     
        m=i;
        
        q=l;
        
        s0=0;
        
        for n1=1:M+1
        
        s0=s0+s(m,n1)*f(u(:,n1,l)); 
        
        end
        
        %s0=s(m,1)*f(u(:,1,l)) + s(m,2)*f(u(:,2,l));
 
% newton loop 
  
           for p=1:max_iter 
        
               Fn1=fun_corr(v,u(:,:,:),dt,m,q,Mass,s0);
        
               %Jn1=jac(v,dt); % analytic jacobian
               
               Jn1=NumJacob(@fun_corr,v,u(:,:,:),dt,m,q,Mass,s0); % numerical jacobian
        
               z=v-Jn1\Fn1; 
        
                   if norm(abs(z-v),inf)<tol
%                
%                    disp('converged')
%                 
                    break
                   end
             v=z;
             
           end
           
       u(:,m+1,l+1) = z;
   
       end
   
    end
   
   y(:,j+1)= u(:,M+1,M+1);  
   
   u0=y(:,j+1);
    
end 

u_st(:,k)=y(:,J+1);

% exact solution at t=T

u_ex=[sin(T);cos(T)];

er(k)=norm( abs(y(:,J+1)-u_ex),inf);

end

p1=log2(er(n-5)/er(n-4));

p2=log2(er(n-4)/er(n-3));

p3=log2(er(n-3)/er(n-2));

p4=log2(er(n-2)/er(n-1));

p5=log2(er(n-1)/er(n));

orders=[0;p1;p2;p3;p4;p5];
   
errors=er;  

% p1=log2(er(n-7)/er(n-6));
% 
% p2=log2(er(n-6)/er(n-5));
% 
% p3=log2(er(n-5)/er(n-4));
% 
% p4=log2(er(n-4)/er(n-3));
% 
% p5=log2(er(n-3)/er(n-2));
% 
% p6=log2(er(n-2)/er(n-1));
% 
% p7=log2(er(n-1)/er(n));
% 
% orders=[0;p1;p2;p3;p4;p5;p6;p7];
%    
% errors=er;


end
