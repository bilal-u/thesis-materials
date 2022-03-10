
function F = fun_F(p,y_old,dt,t)

 F = mass(t,p)*(p-y_old) - dt*f(t,p);
 
end