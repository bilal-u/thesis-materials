
function F = fun_corr(p,t,dt,s1,f )

F= p-dt*f(t,p)- s1;

end
