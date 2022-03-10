
function rh=Rho(y,N)

    u = y(1:N);  
    
    x = y(N+1:end);
    
    rho = zeros(N,1);
    
    v = zeros(N,1);
    
	        for j=2:N-1
                
               v(j)=2/(x(j+1)-x(j-1))*( (u(j+1)-u(j))/(x(j+1)-x(j))-(u(j)-u(j-1))/(x(j)-x(j-1)) );
		   
                end
            
            v(1) = 2*((x(2)-x(1))*(u(3)-u(1))-(x(3)-x(1))*(u(2)-u(1)))/((x(3)-x(1))*(x(2)-x(1))*(x(3)-x(2)));
            
            v(N) = 2*((x(N-1)-x(N))*(u(N-2)-u(N))-(x(N-2)-x(N))*(u(N-1)-u(N)))/((x(N-2)-x(N))*(x(N-1)-x(N))*(x(N-2)-x(N-1)));
            
            rho = rho + v.^2;
            
   % alpha calculation
   
   gamma = 1/3;
   
   Alpha = eps;
   
   for j=2:N
       
      Alpha = Alpha + (1/2)*(rho(j)^gamma+rho(j-1)^gamma)*(x(j)-x(j-1));
      
   end
   
   Alpha = (Alpha)^(3);
   
   % mesh density function
   
   rho = (1+(1/Alpha)*rho).^(1/3);
   
   % smoothing mesh density function
   
   rho_sm=zeros(N,1);
   
   for j=2:(N-1)
       
      rho_sm(j) = 1/4*(rho(j-1)+rho(j+1))+1/2*rho(j);
         
   end
   
   rho_sm(1) = 1/2*(rho(1)+rho(2));
   
   rho_sm(N) = 1/2*(rho(N)+rho(N-1));
   
   rh=rho_sm;
   
  end
