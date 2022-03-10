
function F = fun_corr(p,u,dt,m,l,s0,s1,t )

F=p-u(:,m,l+1)-dt*( inv(mass(t,p))*f(t,p)-s0*f(t,u(:,m+1,l)) )-dt*s1;

end
