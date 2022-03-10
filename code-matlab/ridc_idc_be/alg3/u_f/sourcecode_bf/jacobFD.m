 function J = jacobFD(t)%(g,x,varargin) 
% Calculates the Jacobian of the
% system of non-linear equations:
% f(x) = 0, through finite differences.
% The Jacobian is built by columns
dt=0.01;
J=1.0-dt*( -2.0*pi*sin(2.0*pi*t) - 2.0*( 1.0-cos(2.0*pi*t) ) );

% delx=1e-8;
% 
% m=length(x);
% 
% J=zeros(m,m);
% 
% for j = 1:m
%    xx = x; 
%    xx(j) = x(j) + delx; 
%    
%    f_1=feval(g,x,varargin{:});
%    
%    f_2=feval(g,xx,varargin{:});
%   
%    J(:,j)= (f_2-f_1)/delx; 
end
