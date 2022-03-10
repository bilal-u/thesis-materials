
function rh = integration_matrix(m)

format compact, format long 

syms x;

p=m;

M=p-1;

S=zeros(p-1,p,p-1);

l=1;

while(l<=M)

u=linspace(0,1,l+1);

ln=length(u);

  %for m=1:l
  
  m=1;
  
  while(m<=l)
    

       for i=1:ln
    
       y=1; 
  
       denom=1;
    
       integrant=1;
     
           for k=1:ln
        
              if i~=k 
         
               y=y*(x-u(k));
           
               denom=denom*(u(i)-u(k));
           
              end
         
           end
   
       integrant=integrant*y;
   
       S(m,i,l)=int (integrant,u(m),u(m+1))/denom;

       end
      
   m=m+1;
   
  end

  l=l+1;

end

rh=S;
 
end


% function rh = integration_matrix(m)
% 
% format compact, format long 
% 
% syms x;
% 
% p=m;
% 
% M=p-1;
% 
% S=zeros(p-1,p,p-1);
% 
% l=1;
% 
% while(l<=M)
% 
% u=linspace(0,1,l+1);
% 
% ln=length(u);
% 
%   %for m=1:l
%   
%   m=1;
%   
%   while(m<=l)
%     
% 
%        for i=1:ln
%     
%        y=1; 
%   
%        denom=1;
%     
%        integrant=1;
%      
%            for k=1:ln
%         
%               if i~=k 
%          
%                y=y*(x-u(k));
%            
%                denom=denom*(u(i)-u(k));
%            
%               end
%          
%            end
%    
%        integrant=integrant*y;
%    
%        S(m,i,l)=int (integrant,u(m),u(m+1))/denom;
% 
%        end
%       
%    m=m+1;
%    
%   end
% 
%   l=l+1;
% 
% end
% 
% rh=S;
%  
% end
%  
% 
%  
