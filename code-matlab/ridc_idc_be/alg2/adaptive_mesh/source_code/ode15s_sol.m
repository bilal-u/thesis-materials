
function rh=ode15s_sol(T,N)

h = 1/(N+1);

% initial condition & guess

xint = h*(1:N)';

uint = sin(2*pi*xint) + (1/2)*sin(pi*xint);

u0=[uint; xint];

tspan=[0 T];
  
% option setting  

opts = odeset('RelTol',1e-18,'AbsTol',1e-16,'Mass',@mass_ode,'MaxOrder',5);

sol= ode15s(@f_ode,tspan,u0,opts);

y = deval(sol,T);

rh=y;

end
