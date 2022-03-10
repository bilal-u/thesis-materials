 function J = jacobFD(g,x,varargin) 
% Calculates the Jacobian of the
% system of non-linear equations:
% f(x) = 0, through finite differences.
% The Jacobian is built by columns

delx=1e-8;

m=length(x);

J=zeros(m,m);

for j = 1:m
   xx = x; 
   xx(j) = x(j) + delx; 
   
   f_1=feval(g,x,varargin{:});
   
   f_2=feval(g,xx,varargin{:});
  
   J(:,j)= (f_2-f_1)/delx; 
end
