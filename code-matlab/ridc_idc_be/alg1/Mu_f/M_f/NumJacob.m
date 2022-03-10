function df=NumJacob(f,x0,varargin)

epsilon = 1e-8; % delta

l_x0=length(x0); % length of x0;

f0=feval(f,x0,varargin{:}); % caclulate f0


   for i=1:l_x0
    
      dx = [ zeros(i-1,1); epsilon; zeros(l_x0-i,1)];
    
      df(:,i) = ( feval(f,x0+dx,varargin{:}) - f0)/epsilon;
      
   end
   
end