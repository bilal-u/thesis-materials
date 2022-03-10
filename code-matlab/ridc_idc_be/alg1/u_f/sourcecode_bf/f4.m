
function F=f4(t,u)

   ep=4e-1;

   N=length(u);
   
   h=1/N;

   F=zeros(N,1);
   
   F(1)=0; % boundary condition
      
    for j=2:N-1
           
     F(j)= (ep/h.^2)*( u(j-1) -2*u(j) + u(j+1) );        
          
    end
       
  F(N)=0; % boundary condition
  
end