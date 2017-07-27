%  m = 2*8/4/2*3
%  
 clear global
% % ans = 8/(1/2)*2;
% % guess=[85*pi/180 asin(1/2)+pi/180]
% 
% 
% % 
% % guess=[85*pi/180 asin(1/8.81)+pi/180];    % ????Guess function is used to solve for numerical solution, needto define a interval.
% % beta=fzero(@ThetaBetaSolve,guess,[],6,8.81);
% 
% % eq := ode((1-M^2)/(1+(gamma-1)/(2*y(x)^2)/(gamma*y(x)^2)/(y(x)^2)*y'(x)/(y(x)^2)= 4*f/D, y(x))
% % 
% % solve(eq)
% 
% 
% %x = fzero(x*sin(cosh(x)))
% function y = roughwork()
%     M1 = 3.5;
%     global fm;
%     fm = 0.57643;
%     My  = fsolve(@ffmi,M1)
%     
% end
% 
% function y = ffmi(Mx)
%     global fm;
%     y = -fm + (1-Mx^2)/1.4/Mx^2 + (1.4+1)/2/1.4* log(((1.4+1)*Mx^2)/(2+(1.4-1)*Mx^2));
% end    
%     


%Mout = fanno_flow_2_Modified(3.5,0.5,0.4,2*10^-3)


% x = -pi:pi/20:pi;
% y1 = sin(x);
% plot(x,y1)
% 
% hold on
% y2 = cos(x);
% plot(x,y2)
% hold off
% 
% lgd = legend('sin(x)','cos(x)');
% title(lgd,'My Legend Title')
% x = 0:0.1:10;
% y1 = sin(2*x);
% y2 = cos(2*x);
% 
%          % plot against left y-axis  
% plot(x,y1)           
% hold on       % plot against right y-axis
% plot(x,y2)
% title('Subplot 4')

% x=[2.22594:0.05:2.73436];
%  m=[2.3890 2.3891 2.3782 2.4211 2.6127 2.7343 2.6273 2.4738 2.3638 2.2791 2.2259 2.2516];
%  s=[0.0000 0.0000 0.0912 0.1416 0.1906 0.2200 0.2320 0.2378 0.2414 0.2453 0.2525 0.2791];
%  y=zeros(length(x),length(m));
%  z=zeros(1,12);
%  for i=1:12
%     y(:,i)=normpdf(x,m(1,i),s(1,i));
%     plot(x,y);
%     hold on
%     z(:,i)=prctile(y,[5 25 75 95]);
%  end

gamma = 1.4;
Cf = 0.002;
Dm = 0.06985;
xmax = 2;

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
