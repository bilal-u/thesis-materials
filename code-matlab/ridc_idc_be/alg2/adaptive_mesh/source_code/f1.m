
function rh = f1(t,y)

    N1=length(y);

    N=N1/2;

    ep=1e-2;
    
    tau=1e-2;
    
    u = y(1:N); 
    
    x = y(N+1:end);
    
    x0 = 0;
    
    u0 = 0;
    
    xNP1 = 1;
    
    uNP1 = 0;
    
    g = zeros(2*N,1);
   
    for i = 2:N-1
        
      dx = x(i+1) - x(i-1);
      
      g(i) = (2*ep)/dx*( (u(i+1) - u(i))/(x(i+1) - x(i)) - (u(i) - u(i-1))/(x(i) -...
          x(i-1)) )- (1/2)*(u(i+1)^2 - u(i-1)^2)/dx;
    end
    
    dx = x(2) - x0;    
    
    g(1) = (2*ep)/dx*((u(2) - u(1))/(x(2) - x(1)) - (u(1) - u0)/(x(1) - x0)) -...
        (1/2)*(u(2)^2 - u0^2)/dx;
    
    dx = xNP1 - x(N-1);   
    
    g(N) = (2*ep)/dx*((uNP1 - u(N))/(xNP1 - x(N)) - (u(N) - u(N-1))/(x(N) -...
        x(N-1)))/dx - (1/2)*(uNP1^2 - u(N-1)^2)/dx;
   
    rho_sm=Rho(y,N);

    for i = 2:N-1
        
      g(i+N) = (rho_sm(i+1) + rho_sm(i))*(x(i+1) - x(i)) - ...
               (rho_sm(i) + rho_sm(i-1))*(x(i) - x(i-1));
    end
    
    
    g(1+N) = (rho_sm(2) + rho_sm(1))*(x(2) - x(1)) - (rho_sm(1) + rho_sm(1))*(x(1) - x0);
    
    g(N+N) = (rho_sm(N) + rho_sm(N))*(xNP1 - x(N)) - (rho_sm(N) + rho_sm(N-1))*(x(N) - x(N-1));
        
    g(1+N:end) = - g(1+N:end)/(2*tau);   
    
    rh = g;
    
 end
