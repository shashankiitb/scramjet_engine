clc;
clear;

%variation of properties over combustor length
%only till mach1 choking
%check values of constants in functions

%STATION 4
M4=4;
T04=3676;       %(K)
P4=104*10^3;    %(Pa)

%% TAKE CARE OF VALUES IN THIS SECTION
g=1.4;
R=287;
%COMBUSTOR PROPERTIES
f=0.1;
pbya=4/0.05;  %circular
k=20*10^6;    % heat add rate j/m/kg
%CALCULATING COEFFICIENTS
cp=g*R/(g-1); 
B= f*pbya/2/k;
c1=g;
c2=2*g*cp*B;
c3=(g-1)/2;

%%
T4=T04*(1+(g-1)/2*M4^2)^(-1);
rho4=P4/R/T4;
v4=M4*(g*R*T4)^(0.5);
m_dotbyA=rho4*v4;

%SOLVING ODE
[m52,t05]=ode45('dydxf1',[M4^2,1],T04);
m5=m52.^0.5;

x(:,1)=100*cp*(t05(:,1)-T04)/k;                 %in cms
t5(:,1)=t05(:,1).*(1+(g-1)/2.*m52(:,1)).^(-1);
v5(:,1)=m5(:,1).*(g*R*t5(:,1)).^(0.5);
rho5(:,1)=m_dotbyA./v5(:,1);
p5=rho5.*R.*t5;
p05=p5.*(1+(g-1)/2.*m52(:,1)).^(g/(g-1));

plot(x,m5,'blue');
ylabel('Mach number');
plot(x,t05,'red');
% ylabel('Stagnation Temp');
% plot(x,t5,'green');
% % ylabel('Temperature');
% plot(x,p05,'red');
% ylabel('Stagnation Pressure');
% plot(x,p5,'green');
% ylabel('Pressure');
% plot(x,rho5,'red');
% ylabel('Density');
xlabel('Length');