
function F = fun_corr(p,u,dt,m,l,s0,s1 )

F=p-u(:,m,l+1)-dt*( inv(mass(p))*f(p)-s0*f(u(:,m+1,l)) )-dt*s1;

end