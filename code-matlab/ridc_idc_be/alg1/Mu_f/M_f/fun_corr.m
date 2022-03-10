
function F = fun_corr(p, u,dt,m,l,M,s0)

F= M*(p-u(:,m,l+1))-dt*( f(p)-f(u(:,m+1,l)) ) - dt*s0; 

end