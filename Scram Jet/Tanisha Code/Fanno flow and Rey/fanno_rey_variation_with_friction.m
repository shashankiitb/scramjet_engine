% shows nature of graph_ check if same properties in function dydx
clc;clear;

%STATION 4
M4=4;
T04=3676;

%%Since coefficients not being used in main code here
% g=1.4;
% R=287;
% %COMBUSTOR PROPERTIES
% f=0.001;
% pbya=4/0.05;  %circular
% k=20*10^6;    % heat add rate j/m
% 
% %CALCULATING COEFFICIENTS
% cp=g*R/(g-1); 
% B= f*pbya/2/k;
% c1=g;
% c2=2*g*cp*B;
% c3=(g-1)/2;
gamma = 1.4;
Cf = 0.002;
Dm = 2.75;
xmax = 22.78;
f = @(t,x)[-(x(1)/2)*(1+((gamma-1)*(x(1)*x(1)))/2)*(93*Cf/(Dm*x(2))); Cf/(2*Dm)*(93+x(1)*x(1)*(93*(gamma-1)-89*gamma*x(2))); 89/Dm*Cf*gamma*(x(1)*x(1))*x(3)/2];

%SOLVING ODE
[t1,xa1] = ode45(f,[0,xmax],[2.60,1.0,2.61]);
plot(t1,xa1(:,3),'red');
hold on;
[t2,xa2] = ode45(f,[0,xmax],[2.60,1.0,3.103]);
plot(t2,xa2(:,3),'blue');
hold on;
[t3,xa3] = ode45(f,[0,xmax],[2.60,1.0,3.509]);
plot(t3,xa3(:,3),'green');
hold on;

% %% REYLEIGH FLOW
% %Reyleigh differential solver
% [m42_rey,t04_rey]=ode45('dydx_rey',[M4^2,0.625],T04);
% scatter(m42_rey.^0.5,t04_rey,'black');
% hold on;
% 
% legend('f=0.1','f=0.01','f=0.001','f=0');
%Reyleigh actual non diff equation
% t_rey2(:,1)=T04*Rey(m42_rey(:,1).^(0.5),g)/Rey(M4,g);
% scatter(m42_rey.^0.5,t_rey2,'blue');


% %To check for last point, which arent plotted for some reason
% scatter(m42_f1(41,1),t04_f1(41,1));
% scatter(m42_f2(45,1),t04_f2(45,1));
% scatter(m42_f3(45,1),t04_f3(45,1));