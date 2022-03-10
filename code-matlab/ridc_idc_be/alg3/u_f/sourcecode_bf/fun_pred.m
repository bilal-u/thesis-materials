
function F = fun_pred(p,u_old,t,dt,f) 

F= p-dt*f(p)-u_old; 

end