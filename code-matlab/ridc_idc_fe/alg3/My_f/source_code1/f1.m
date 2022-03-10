
function rh=f1(t,x)

% Mass matrix 

L=[4 -1;
    -1 4];
% rhs 

f=[x(1)+4*x(2);
    -4*x(1)-x(2)];

% linear solve

h=gauss(L,f);

%h=L\f;

rh= h;

end
