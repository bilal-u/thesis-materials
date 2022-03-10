function rh=beuler(y0,tspan,dt)

N=length(y0);

T=tspan(1,2);

W=int16(T/dt);

u=zeros(N,W+1);

u(:,1)=y0;

max_iter=100;

tol=1e-8;

it_count=zeros(W,1);

  for m=1:W
        
        v=u(:,m);
     
        it=0;

      for p=1:max_iter 
        
        Fn=fun_F(v,u(:,m),dt);
        
        Jn=jacobFD(@fun_F,v,u(:,m),dt); 
       
      % Jn=NumJacob(@fun_F,v,u(:,m-1),dt);
        
        z=v-Jn\Fn; 
        
        it=it+1;
       
           if norm(abs(z-v))<tol
          
           %  disp('converged')
             break
             
           end
           
        v=z;
        
      end
    
    u(:,m+1) = z;
    
    it_count(m)=it;
    
  end    
   
rh=u(:,end);

%it_count

end