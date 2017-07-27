clc;clear;
gamma = 1.4; %%% ???
Cf = 0.002;  %%% ???
Dm = 0.012;  %%% ??? in meters
xmax = 2;    

%%%%%%%%%%%%%%%%%%%%%%%%

% This formula is taken from, assuming epslion (roughness factor) = 0
% https://en.wikipedia.org/wiki/Fanning_friction_factor 
% Re = 
% A = (2.457*log((Re/7)^0.9))^16;
% B = (37530/Re)^16;
% Cf = 2((8/Re)^12 + (A+B)^(-1.5))^(1/12);


%%%%%%%%%%%%%%%%%%%%%%%%

f = @(t,x)[-(x(1)/2)*(1+((gamma-1)*(x(1)*x(1)))/2)*(93*Cf/(Dm*x(2))); Cf/(2*Dm)*(93+x(1)*x(1)*(93*(gamma-1)-89*gamma*x(2))); 89/Dm*Cf*gamma*(x(1)*x(1))*x(3)/2];

[t1,xa1] = ode45(f,[0,xmax],[2.231,1.0,36515]);

plot(t1,xa1(:,3));
xlabel('x(m)'), ylabel('Pw');
grid on
grid minor
figure;
plot(t1,xa1(:,1));
xlabel('x(m)'), ylabel('M');
grid on
grid minor


lgdx = legend(' 45.0 ');