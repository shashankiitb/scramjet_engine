clc;clear;
gamma = 1.4;
Cf = 0.002;
Dm = 2.75;    %inches
xmax = 22.78-10.49; %inches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%syms M(x) R(x) P(x)

% Defining Diff Equation   
% x(1) = M(x), x(2) = R(x), x(3) = P(x)
%ode1 = diff(x(1)) == -(x(1)/2)*(1+((gamma-1)*(x(1)*x(1)))/2)*(93*Cf/(Dm*x(2)));
%ode2 = diff(x(2)) == Cf/(2*Dm)*(93+x(1)*x(1)*(93*(gamma-1)-89*gamma*x(2)));
%ode3 = diff(x(3)) == 89/Dm*Cf*gamma*(x(1)*x(1))*x(3)/2;
%odes = [-(x(1)/2)*(1+((gamma-1)*(x(1)*x(1)))/2)*(93*Cf/(Dm*x(2))); Cf/(2*Dm)*(93+x(1)*x(1)*(93*(gamma-1)-89*gamma*x(2))); 89/Dm*Cf*gamma*(x(1)*x(1))*x(3)/2]

f = @(t,x)[-(x(1)/2)*(1+((gamma-1)*(x(1)*x(1)))/2)*(93*Cf/(Dm*x(2))); Cf/(2*Dm)*(93+x(1)*x(1)*(93*(gamma-1)-89*gamma*x(2))); 89/Dm*Cf*gamma*(x(1)*x(1))*x(3)/2];
%Above is main function which consists all of the main 3 odes as mention in above comment 

%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting graph for different Mach No.
%%%%%%%%%%%%%%%%%%%%%%%%
mstart = 2.60;    % Starting point for Mach no.
mend = 4;         % End point for mach no.
mstep = 0.1;      % difference between subsequent Mach no.
Min = mstart:mstep:mend;    % Mach no. Array
iter = round(((mend-mstart)/mstep) +1)  % No of iteration or size of mach no array
ti = zeros(45, iter);
xa = zeros(45,3,iter);

 for i=1:iter
     [ti(1:45,i),xa(:,:,i)] = ode45(f,[0,xmax],[Min(1,i),1.0,2.3625]);
     % here initial condition is that R(x) is 1, P is 2.3625(inital
     % pressure value obtained from bilig's paper)
 end

a1 = (xa(:,3,:))/45.0; % 45.0 Psia is a inital total pressure Pt0 


for i=1:iter
    plot(ti(:,i),a1(:,:,i));
    hold on;
end

grid on
grid minor
xlabel('x(Inches)'), ylabel('Pw/Pt0');
lgdx = legend('2.60 = M','2.7 = M','2.80 = M','2.90 = M','3.0 = M','...');