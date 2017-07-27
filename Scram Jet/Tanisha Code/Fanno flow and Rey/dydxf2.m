function dydxf2=dydxf2(x,y)   %x=M^2, y=T0

%match these values as in code
g=1.4;
R=287;
%COMBUSTOR PROPERTIES
f=0.01;
pbya=4/0.05;  %circular
k=20*10^6;    % heat add rate j/m

%CALCULATING COEFFICIENTS
cp=g*R/(g-1); 
B= f*pbya/2/k;
c1=g;
c2=2*g*cp*B;
c3=(g-1)/2;

dydxf2=y/x*(1-x)/(1+c3*x)/(1+c1*x+c2*x*y);
end