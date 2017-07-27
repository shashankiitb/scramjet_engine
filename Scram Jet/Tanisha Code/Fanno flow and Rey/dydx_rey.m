%reyleigh function solver
function dydx_rey=dydx_rey(x,y)   %x=M^2, y=T0
g=1.4;
c1=g;
c3=(g-1)/2;
dydx_rey=y/x*(1-x)/(1+c3*x)/(1+c1*x);
end