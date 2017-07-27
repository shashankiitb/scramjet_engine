% Below Code give us scramjet analysis from 


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

for e = 1:iter
    %Analysis for Beta solution.
    for i = 1:3
        guess=[asin(1/M(i,e))+10^-6];    % Giving inital point for fzero, as (Mnormal should be > 1)
        beta(i,e)=fzero(@ThetaBetaSolve,guess,[],theta(i,e),M(i,e));  %using function fzero to solve nonlinear equation
      
%using oblique shock relations to find downstream properties
        M(i+1,e) = (((gamma_a-1)*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))+2)/(2*gamma_a*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))-(gamma_a-1))/(sin(beta(i,e)-theta(i,e))*sin(beta(i,e)-theta(i,e))))^(1/2);
        p(i+1,e) = p(i,e)*(2*gamma_a*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))-(gamma_a-1))/(gamma_a+1);
        rho(i+1,e) = rho(i,e)*((gamma_a+1)*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e)))/((gamma_a-1)*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))+2);
        t(i+1,e) = p(i+1,e)/p(i,e)*rho(i,e)/rho(i+1,e)*t(i,e);
        v(i,e) = M(i,e)*(gamma_a*R_a*t(i,e))^(1/2);
        t0(i+1,e) = t(i+1,e)*(1+(gamma_a-1)/2*M(i+1,e)*M(i+1,e));
        p0(i+1,e) = p(i+1,e)*(t0(i+1,e)/t(i+1,e))^(gamma_a/(gamma_a-1));
    end

   %Calculation of inlet parameters
    inlet_drag(1,e)=(p(2,e)*h1*w)+(p(3,e)*h2*w);
    inlet_lift(1,e)=(p(2,e)*h1*w*cot(theta(1,e)))+(p(3,e)*h2*w*cot(theta(3,e)))-(p(3,e)*h3*cot(beta(3,e)-theta(3,e))*w);
    totalpressureloss(1,e)=p0(4,e)/p0(1,e); 

    %Combustor Properties
    MWair = 28.85; % Molecular weight of air
    MWfuel = 2.00;  % Molecular Weight of fuel ( here H2 )
    a= 0.5; % For H2 as comustor material
    nstoic = 4.76 * a;
    mstoic = nstoic * MWair/MWfuel;
    effcom = 1; % Effecincy of combustor
    mdot(1,e)=rho(1,e)*v(1,e)*w*hc;
    mdot_fuel(1,e)=phi*mdot(1,e)/mstoic; %using definition of equivalence ratio to calculate fuel massflow rate, 0.06796 is the stoichiometric air-fuel mass flow ratio.
    lhv=120*10^6;  %lower heating value of fuel
    q(1,e)=mdot_fuel(1,e)*lhv*effcom/mdot(1,e); %specific heat to be added

    %Rayleigh flow
    t0(5,e)=t0(4,e) + q(1,e)/cpa;

    t0max=t0(4,e)/(2*(1+gamma_a)*M(4,e)*M(4,e))/(1+(gamma_a-1)/2*M(4,e)*M(4,e))*(1+gamma_a*M(4,e)*M(4,e))^(2); %calculation of critical total temperature
    
%using bisection method to solve for M5 ????????
    a=1;
    b=M(4,e);
    M(5,e)=(a+b)/2;
    for j=1:50
        LHS=t0max/t0(5,e);
        RHS=1/(2*(1+gamma_g)*M(5,e)*M(5,e))/(1+(gamma_g-1)/2*M(5,e)*M(5,e))*(1+gamma_g*M(5,e)*M(5,e))^(2);
        if(LHS-RHS)>0
            a=M(5,e);
            M(5,e)=(a+b)/2;
        end
        if(LHS-RHS)<0
            b=M(5,e);
            M(5,e)=(a+b)/2;
        end
        if LHS==RHS
            M(5,e)=(a+b)/2;
            break;
        end
    end
    
%calculation of properties in region 5
    p(5,e)=p(4,e)*(1+gamma_a*M(4,e)*M(4,e))/(1+gamma_a*M(5,e)*M(5,e));
    t(5,e)=t0(5,e)/(1+(gamma_a-1)/2*M(5,e)*M(5,e));
    p0(5,e)=p(5,e)*(1+(gamma_a-1)/2*M(5,e)*M(5,e))^(gamma_a/(gamma_a-1));


    %Nozzle Analysis
    A_comb=h3*w;
    %choking area
    A_crit=A_comb/(1/M(5,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(5,e)*M(5,e)))^((gamma_a+1)/2/(gamma_a-1))); 
    A_exit=he*w;

    guess_Noz=[1 30];

    M(6,e)=fzero(@Nozzle,guess_Noz,[],A_exit,A_crit);  %using function fzero to solve for nonlinear equation
    t(6,e)=t(5,e)*(cpa + gamma_a*M(5,e)*M(5,e)*R_a/2)/(cpa + gamma_a*M(6,e)*M(6,e)*R_a/2);
    p(6,e)=p(5,e)*(t(6,e)/t(5,e))^(gamma_a/(gamma_a-1));
    p0(6,e)=p(6,e)*(1+(gamma_a-1)/2*M(6,e)*M(6,e))^(gamma_a/(gamma_a-1));
    t0(6,e)=t0(5,e);

    v(6,e)=M(6,e)*(gamma_a*R_a*t(6,e))^(1/2); 
    v(5,e)=M(5,e)*(gamma_a*R_a*t(5,e))^(1/2);

%calculation of critical performance parameters
    net_thrust(1,e)=(mdot(1,e))*v(6,e) - mdot(1,e)*v(1,e)+p(6,e)*he*w-p(1,e)*hc*w;
    nozzle_thrust(1,e)=(mdot(1,e)+mdot_fuel(1,e))*(v(6,e) - v(5,e))+p(6,e)*he*w-p(5,e)*h3*w;%Pressure Thrust
    g0=9.812;
    isp(1,e)=(net_thrust(1,e))/g0/mdot_fuel(1,e);
    
%calculation of efficiency parameters
    t_eff(1,e)=((mdot(1,e)+mdot_fuel(1,e))*v(6,e)*v(6,e)/2 - mdot(1,e)*v(1,e)*v(1,e)/2)/mdot_fuel(1,e)/lhv; %thermal efficiency
    p_eff(1,e)=2*v(1,e)/v(6,e)/(1+v(1,e)/v(6,e)); %propulsive efficiency
    o_eff(1,e)=t_eff(1,e)*p_eff(1,e); %overall efficiency
end

