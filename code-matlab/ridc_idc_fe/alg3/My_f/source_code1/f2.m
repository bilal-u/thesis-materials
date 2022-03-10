function rh=f2(t,x)

% Mass matrix 

L=[x(1)^2+4 -1/2;
    -1/2 x(2)^2+4];

f=[x(1)^2*x(2)+4*x(2)+(1/2)*x(1);
    -(1/2)*x(2)-x(1)*x(2)^2-4*x(1)];

% linear solve

h=gauss(L,f);

%h=L\f;

rh= h;

end
