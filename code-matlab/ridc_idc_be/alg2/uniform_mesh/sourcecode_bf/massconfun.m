
function F=f3(t,x)

m=[4 -1; 
  -1 4];

g=[x(1)+4*x(2) ;
    -4*x(1)-x(2)];

F=m\g;

end
