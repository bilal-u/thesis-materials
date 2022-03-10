
function rh=f(t,y)

 N=length(y);
 
 ep=1e-2;
 
 h=1/(N-1);
 
 yp=zeros(N,1);
 
 yp(1)=0; % BC
 
 for j=2:N-1
     
 yp(j)=(ep/h^2)* (y(j+1)-2*y(j)+y(j-1)) - 1/(4*h)* (y(j+1)^2-y(j-1)^2) ;
 
 end
 
 yp(N)=0; %Bc
 
 rh=yp;

end
