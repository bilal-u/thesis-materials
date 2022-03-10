
function F = fun_pred(u,u_old,dt,M) 

F= M*(u-u_old)-dt*f(u); 

end