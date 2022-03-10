
function rh = integration_matrix(m)

format compact, format long 

syms x;

p=m;

M=p-1;

u=linspace(0,1,p);

l=length(u);

m=1;

S=zeros(M,p);

while(m<=M)

for i=1:l
    
   % m=1;
    
    y=1; 
  
    denom=1;
    
    integrant=1;
     
    for k=1:l
        
        if i~=k 
         
           y=y*(x-u(k));
           
           denom=denom*(u(i)-u(k));
           
        end
        
    end
 
   
integrant=integrant*y;
   
S(m,i)=int (integrant,u(m),u(m+1))/denom;
   
    
end

 m=m+1;
 
 rh=S;
 
 end
