clear global

iter = 1; % No.of Iteration of Mach No.
%initial conditions
Mdes=7;         %design point Mach number
alt=26;            %altitude: 26 km 
R_a=287.1;         %Specific gas constant of air
R_g=287.1;         %Specific gas constant of air-fuel ratio
p1=2188;           %in Pa
t1=222.5;          %in K
rho1=p1/R_a/t1;    %P=rho*R*T
gamma_a=1.4;       %ratio of specific heats, for air
gamma_g=1.4;       %ratio of specific heats, for air-fuel mixture
cpa = 1003; %gamma_a*R_a/(gamma_a-1);     %gamma*R/(gamma-)


%body parametersnn
h1=0.3364;
h2=0.5794;
h3=0.0842;
hc=h1+h2+h3;    %capture height
he=1.5;         %exit height
w=1;            %width of engine

beta=zeros(3,iter);     %array for shock angle    
theta=zeros(3,iter);    %array for turn angle, which is same for all freestream Mach numbers

phi=1;         %equivalence ratio of engine

%defining fluid property matrices
M=zeros(6,iter);       %Mach number
v=zeros(6,iter);       %Velocity
p=zeros(6,iter);       %Pressure
t=zeros(6,iter);       %Temperature
p0=zeros(6,iter);      %Total Pressure
t0=zeros(6,iter);      %Total Temperature
rho=zeros(6,iter);     %Density

%The first index signifies the subdivision number in the engine, while the second 
%defining arrrays for different parameters ??? Need to add more comment

inlet_drag=zeros(1,iter);      %inlet drag
inlet_lift=zeros(1,iter);      %inlet lift
totalpressureloss=zeros(1,iter);    %total pressure loss in inlet, due to shocks

nozzle_thrust=zeros(1,iter);   %Thrust produced by only the nozzle
net_thrust=zeros(1,iter);      %Net force, or installed thrust
isp=zeros(1,iter);             %Specific Impulse

%defining shock angles
for i = 1:iter
    theta(1,i)=6/180*pi;
    theta(2,i)=19.705/180*pi;
    theta(3,i)=25.705/180*pi;   
    
    M(1,i)=Mdes;%*(100+i-11)/100;  %definition of Mach number range
    
%defining initial conditions
    p(1,i)=p1;
    t(1,i)=t1;
    rho(1,i)=rho1;
    v(1,i)=M(1,i)*(gamma_a*R_a*t(1,i))^(1/2);
    t0(1,i)=t(1,i)*(1+(gamma_a-1)/2*M(1,i)*M(1,i));
    p0(1,i)=p(1,i)*(t0(1,i)/t(1,i))^(gamma_a/(gamma_a-1));
end
