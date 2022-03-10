
function [orders,errors]=idc1(a,T)

M=1;

dt_ar=[a;a/2;a/4;a/8;a/16;a/32];

n=length(dt_ar);

% error storing

er=zeros(n,1)';

for k=1:n
    
dt=dt_ar(k);

W=T/dt; % number of subinterval for t

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
    
 u_st(:,j+1)= Y(:,M+1,1);  
   
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