
function F = fun_F(p,y_old,dt,t,f)

 F = mass(t,p)*(p-y_old) - dt*f(t,p);
 
end