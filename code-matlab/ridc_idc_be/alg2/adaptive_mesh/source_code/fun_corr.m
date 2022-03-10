
function F = fun_corr(p,u,dt,t,m,l,s0,s1,f )

% v2,Y(:,:,:),dt,t(m+1),i,q,s2,s3,f
F=p-u(:,m,l+1)-dt*( inv(mass(t,p))*f(t,p)-s0 )-dt*s1; 

end