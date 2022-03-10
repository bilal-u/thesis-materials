
function F = fun_F(p,y_old,dt)

 F = mass(p)*(p-y_old) - dt*f(p);
 
end