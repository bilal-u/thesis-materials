
function F = fun_pred(p,u_old,t,dt,f) 

F= p-u_old-dt*f(t,p); 

end