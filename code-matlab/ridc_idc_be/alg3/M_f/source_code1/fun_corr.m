
function F = fun_corr(p,u,t,dt,s1,m,l,f,M )

F= p-dt*( f(t,p) -f(t,u(:,m+1,l)) ) - u(:,m,l+1) - M*dt*s1; 

end
