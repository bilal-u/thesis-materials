
function F=f6(t,u)

  ep=1e-2;   

  N=length(u);
   
  h=1/N;

  F=zeros(N,1);
         
    for j=2:N-1
           
     F(j)= (ep/h^2)*( u(j-1) -2*u(j) + u(j+1) ) + 1/(2*h)*( u(j+1)-u(j-1) );        
          
    end
       
  F(1)=(ep/h^2)*( u(N-1) -2*u(1) + u(2) ) + 1/(2*h)*( u(2)-u(N-1) ); % BC periodic u0=u
  
  F(N)=(ep/h^2)*( u(N-1) -2*u(N) + u(2) ) + 1/(2*h)*( u(2)-u(N-1) );
  
end