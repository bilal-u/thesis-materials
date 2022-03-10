 
function rh = mass(t,y) 

    N1=length(y);

    N=N1/2;

    u = y(1:N);  
    
    x = y(N+1:end);
    
    x0 = 0;
    
    u0 = 0;
    
    xNP1 = 1;
    
    uNP1 = 0;
    
    M1 = speye(N);
    
    M2 = sparse(N,N);
    
    M2(1,1) = - (u(2) - u0)/(x(2) - x0);
    
         for i = 2:N-1
             
            M2(i,i) = - (u(i+1) - u(i-1))/(x(i+1) - x(i-1));
            
         end
         
    M2(N,N) = - (uNP1 - u(N-1))/(xNP1 - x(N-1));
    
    M3 = sparse(N,N);
    
    e = ones(N,1);
    
    M4 = spdiags([e -2*e e],-1:1,N,N);
    
    rh = [M1 M2
          M3 M4];
      
  end
 