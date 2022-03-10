

function out = f(y)
  
    N1 = length(y);
   
    N=N1/2;
   
    ep = 5e-1;
    
    tau = 1e-1;
  
    u = y(1:N);  
    
    x = y(N+1:end);
    
    L = 10;
    
    x_0 = L/2;

    sg = 1/2;
    
    x0 = 0;
    
    u0 = exp(-((0-x_0)^2)/sg); 
    
    xNP1 = 10;
    
    uNP1 =  exp(-((10-x_0)^2)/sg);
    
    g = zeros(2*N,1);
    
    for i = 2:N-1
       g(i) = ep*((u(i+1) - u(i-1))/(x(i+1) - x(i-1))); 
    end
    
    g(1) = ep*((u(2) - u0)/(x(2) - x0) );      
        
    g(N) = ep*((uNP1 - u(N))/(xNP1 - x(N)) );

    % Evaluate monitor function.
    
    M = zeros(N,1);
    
    for i = 2:N-1
       M(i) = sqrt(1 + ((u(i+1) - u(i-1))/(x(i+1) - x(i-1)))^2);
    end
    
    M0 = sqrt(1 + ((u(1) - u0)/(x(1) - x0))^2);
     
    M(1) = sqrt(1 + ((u(2) - u0)/(x(2) - x0))^2);
  
    M(N) = sqrt(1 + ((uNP1 - u(N-1))/(xNP1 - x(N-1)))^2);
    
    MNP1 = sqrt(1 + ((uNP1 - u(N))/(xNP1 - x(N)))^2);
    
    % Spatial smoothing with gamma = 2, p = 2.
    
    SM = zeros(N,1);
    
    for i = 3:N-2
        
      SM(i) = sqrt((4*M(i-2)^2 + 6*M(i-1)^2 + 9*M(i)^2 + ...
                    6*M(i+1)^2 + 4*M(i+2)^2)/29);
    end
    
    SM0 = sqrt((9*M0^2 + 6*M(1)^2 + 4*M(2)^2)/19);
    
    SM(1) = sqrt((6*M0^2 + 9*M(1)^2 + 6*M(2)^2 + 4*M(3)^2)/25);
    
    SM(2) = sqrt((4*M0^2 + 6*M(1)^2 + 9*M(2)^2 + 6*M(3)^2 + 4*M(4)^2)/29);
    
    SM(N-1) = sqrt((4*M(N-3)^2 + 6*M(N-2)^2 + 9*M(N-1)^2 + 6*M(N)^2 + 4*MNP1^2)/29);
    
    SM(N) = sqrt((4*M(N-2)^2 + 6*M(N-1)^2 + 9*M(N)^2 + 6*MNP1^2)/25);
    
    SMNP1 = sqrt((4*M(N-1)^2 + 6*M(N)^2 + 9*MNP1^2)/19);
    
    for i = 2:N-1
      g(i+N) = ((SM(i+1) + SM(i))*(x(i+1) - x(i)) - (SM(i) + SM(i-1))*(x(i) - x(i-1)));
    end
    
    g(1+N) = ((SM(2) + SM(1))*(x(2) - x(1)) - (SM(1) + SM0)*(x(1) - x0));
    
    g(N+N) = ((SMNP1 + SM(N))*(xNP1 - x(N)) - (SM(N) + SM(N-1))*(x(N) - x(N-1)));
    
    g(1+N:end) = - g(1+N:end)/(2*tau);
   
    out = g;
    
  end  

