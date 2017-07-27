%initial conditions
Mdes=9.52; %this is the design mach number, all quantities will be tabulated against 
%Mach number varying by +- of 10% around this design point Mach number.

alt=26; %altitude: 26 km 
R_a=287.1;%specific gas constant of air
R_g=284.09; %specific gas constat for air-fuel mixture in combustor
p1=2163.2;  %ambient pressure, in Pa
t1=219.34; %ambient temperature, in K
rho1=p1/R_a/t1; %ambient air density, in kg/m^3
gamma_a=1.4; %ratio of specific heats for air
gamma_g=1.33; %ratio of specific heats for air-fuel mixture in combustor
cp=gamma_a*R_a/(gamma_a-1); %specific heat of air at constant pressure

%body parameters
h1=0.2855; %height of first ramp of inlet
h2=0.63; %height of second ramp of inlet
h3=0.0845; %combustor height
hc=h1+h2+h3; %capture height
he=1.5;      %exit height
w=1;         %width of engine

beta=zeros(3,21); %wave angle, for each freestream Mach Number
theta=zeros(3,21); %turn angle, for each freestream Mach Number

phi=0.3; %equivalence ratio

%defining fluid property matrices
M=zeros(8,21); %Mach number of flow at 8 points in the engine, for each freestream Mach number
v=zeros(8,21); %Flow velocity at 8 points in the engine, for each freestream Mach number
p=zeros(8,21); %Flow pressure at 8 points in the engine, for each freestream Mach number
t=zeros(8,21); %Flow temperature at 8 points in the engine, for each freestream Mach number
p0=zeros(8,21); %Flow total pressure of flow at 8 points in the engine, for each freestream Mach number
t0=zeros(8,21);%Flow total temperature at 8 points in the engine, for each freestream Mach number
rho=zeros(8,21); %Flow density at 8 points in the engine, for each freestream Mach number

%defining arrrays for different parameters
inlet_drag=zeros(1,21); %pressure drag in inlet, for each freestream Mach number
inlet_lift=zeros(1,21); %pressure lift in inlet, for each freestream Mach number
totalpressureloss=zeros(1,21); %total pressure loss in inlet, due to oblique shocks
A_crit=zeros(3,21); %Critical area in all places of isentropic flow
nozzle_thrust=zeros(1,21); %Thrust produced by nozzle, for each freestream Mach number
net_thrust=zeros(1,21); %Net thrust produced 
isp=zeros(1,21); %Specific impulse of engine, for reach freestream Mach number

t_eff=zeros(1,21); %Thermal efficiency
p_eff=zeros(1,21); %Propulsive efficiency
o_eff=zeros(1,21); %Overall efficiency

for i = 1:21
    %turn angle is same since inlet geometry is fixed
    theta(1,i)=4/180*pi;
    theta(2,i)=8/180*pi;
    theta(3,i)=12/180*pi;   
    
    %establishing the freestream Mach number range, and intial conditions
    %for all freestream Mach numbers
    M(1,i)=Mdes*(100+i-11)/100;
    p(1,i)=p1;
    t(1,i)=t1;
    rho(1,i)=rho1;
    v(1,i)=M(1,i)*(gamma_a*R_a*t(1,i))^(1/2);
    t0(1,i)=t(1,i)*(1+0.2*M(1,i)*M(1,i));
    p0(1,i)=p(1,i)*(t0(1,i)/t(1,i))^(1.4/0.4);
end

for e = 1:21
    %Analysis
    
    for i = 1:3
        guess=[85*pi/180 asin(1/M(i,e))+pi/180]; %initial guess
    
        beta(i,e)=fzero(@ThetaBetaSolve,guess,[],theta(i,e),M(i,e)); %calling function with the expression
    
        M(i+1,e)=(((gamma_a-1)*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))+2)/(2*gamma_a*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))-(gamma_a-1))/(sin(beta(i,e)-theta(i,e))*sin(beta(i,e)-theta(i,e))))^(1/2);
        p(i+1,e)=p(i,e)*(2*gamma_a*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))-(gamma_a-1))/(gamma_a+1);
        rho(i+1,e)=rho(i,e)*((gamma_a+1)*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e)))/((gamma_a-1)*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))+2);
        t(i+1,e)=p(i+1,e)/p(i,e)*rho(i,e)/rho(i+1,e)*t(i,e);
        v(i,e)=M(i,e)*(gamma_a*R_a*t(i,e))^(1/2);
        t0(i+1,e)=t(i+1,e)*(1+0.2*M(i+1,e)*M(i+1,e));
        p0(i+1,e)=p(i+1,e)*(t0(i+1,e)/t(i+1,e))^(1.4/0.4);
    end

    inlet_drag(1,e)=(p(2,e)*h1*w)+(p(3,e)*h2*w);
    inlet_lift(1,e)=(p(2,e)*h1*w*cot(theta(1,e)))+(p(3,e)*h2*w*cot(theta(3,e)))-(p(3,e)*h3*cot(beta(3,e)-theta(3,e))*w);
    totalpressureloss(1,e)=p0(4,e)/p0(1,e);

    %Combustor Properties
    mdot(1,e)=rho(1,e)*v(1,e)*w*hc;
    mdot_fuel(1,e)=phi*0.0679615815*mdot(1,e);
    lhv=120*10^6;
    q(1,e)=mdot_fuel(1,e)*lhv/mdot(1,e);
    
    %Diffuser 1
    
    A_comb_in=h3*w;
    A_crit(1,e)=A_comb_in/(1/M(4,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(4,e)*M(4,e)))^((gamma_a+1)/2/(gamma_a-1)));
    
    M(6,e)=1;

    t0(5,e)=t0(4,e);
t0(6,e)=t0(5,e)* + q(1,e)/(gamma_a*R_g/(gamma_a-1));
    
    a=0.01;
    b=1;
    M(5,e)=(a+b)/2;
    for j=1:50
        LHS=t0(6,e)/t0(5,e);
        RHS=(1+(gamma_a-1)/2*M(6,e)*M(6,e))/(1+(gamma_a-1)/2*M(5,e)*M(5,e))*(M(6,e)*(1+gamma_a*M(5,e)*M(5,e))/M(5,e)/(1+gamma_a*M(6,e)*M(6,e)))^(2);
        if(LHS-RHS)<0
            a=M(5,e);
            M(5,e)=(a+b)/2;
        end
        if(LHS-RHS)>0
            b=M(5,e);
            M(5,e)=(a+b)/2;
     
%defining fluid property matrices
M=zeros(8,21); %Mach number of flow at 8 points in the engine, for each freestream Mach number
v=zeros(8,21); %Flow velocity at 8 points in the engine, for each freestream Mach number
p=zeros(8,21); %Flow pressure at 8 points in the engine, for each freestream Mach number
t=zeros(8,21); %Flow temperature at 8 points in the engine, for each freestream Mach number
p0=zeros(8,21); %Flow total pressure of flow at 8 points in the engine, for each freestream Mach number
t0=zeros(8,21);%Flow total temperature at 8 points in the engine, for each freestream Mach number
rho=zeros(8,21); %Flow density at 8 points in the engine, for each freestream Mach number

%defining arrrays for different parameters
inlet_drag=zeros(1,21); %pressure drag in inlet, for each freestream Mach number
inlet_lift=zeros(1,21); %pressure lift in inlet, for each freestream Mach number
totalpressureloss=zeros(1,21); %total pressure loss in inlet, due to oblique shocks
A_crit=zeros(3,21); %Critical area in all places of isentropic flow
nozzle_thrust=zeros(1,21); %Thrust produced by nozzle, for each freestream Mach number
net_thrust=zeros(1,21); %Net thrust produced 
isp=zeros(1,21); %Specific impulse of engine, for reach freestream Mach number

t_eff=zeros(1,21); %Thermal efficiency
p_eff=zeros(1,21); %Propulsive efficiency
o_eff=zeros(1,21); %Overall efficiency

for i = 1:21
    %turn angle is same since inlet geometry is fixed
    theta(1,i)=4/180*pi;
    theta(2,i)=8/180*pi;
    theta(3,i)=12/180*pi;   
    
    %establishing the freestream Mach number range, and intial conditions
    %for all freestream Mach numbers
    M(1,i)=Mdes*(100+i-11)/100;
    p(1,i)=p1;
    t(1,i)=t1;
    rho(1,i)=rho1;
    v(1,i)=M(1,i)*(gamma_a*R_a*t(1,i))^(1/2);
    t0(1,i)=t(1,i)*(1+0.2*M(1,i)*M(1,i));
    p0(1,i)=p(1,i)*(t0(1,i)/t(1,i))^(1.4/0.4);
end

for e = 1:21
    %Analysis
    
    for i = 1:3
        guess=[85*pi/180 asin(1/M(i,e))+pi/180]; %initial guess
    
        beta(i,e)=fzero(@ThetaBetaSolve,guess,[],theta(i,e),M(i,e)); %calling function with the expression
    
        M(i+1,e)=(((gamma_a-1)*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))+2)/(2*gamma_a*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))-(gamma_a-1))/(sin(beta(i,e)-theta(i,e))*sin(beta(i,e)-theta(i,e))))^(1/2);
        p(i+1,e)=p(i,e)*(2*gamma_a*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))-(gamma_a-1))/(gamma_a+1);
        rho(i+1,e)=rho(i,e)*((gamma_a+1)*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e)))/((gamma_a-1)*M(i,e)*M(i,e)*sin(beta(i,e))*sin(beta(i,e))+2);
        t(i+1,e)=p(i+1,e)/p(i,e)*rho(i,e)/rho(i+1,e)*t(i,e);
        v(i,e)=M(i,e)*(gamma_a*R_a*t(i,e))^(1/2);
        t0(i+1,e)=t(i+1,e)*(1+0.2*M(i+1,e)*M(i+1,e));
        p0(i+1,e)=p(i+1,e)*(t0(i+1,e)/t(i+1,e))^(1.4/0.4);
    end

    inlet_drag(1,e)=(p(2,e)*h1*w)+(p(3,e)*h2*w);
    inlet_lift(1,e)=(p(2,e)*h1*w*cot(theta(1,e)))+(p(3,e)*h2*w*cot(theta(3,e)))-(p(3,e)*h3*cot(beta(3,e)-theta(3,e))*w);
    totalpressureloss(1,e)=p0(4,e)/p0(1,e);

    %Combustor Properties
    mdot(1,e)=rho(1,e)*v(1,e)*w*hc;
    mdot_fuel(1,e)=phi*0.0679615815*mdot(1,e);
    lhv=120*10^6;
    q(1,e)=mdot_fuel(1,e)*lhv/mdot(1,e);
    
    %Diffuser 1
    
    A_comb_in=h3*w;
    A_crit(1,e)=A_comb_in/(1/M(4,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(4,e)*M(4,e)))^((gamma_a+1)/2/(gamma_a-1)));
    
    M(6,e)=1;

    t0(5,e)=t0(4,e);
t0(6,e)=t0(5,e)* + q(1,e)/(gamma_a*R_g/(gamma_a-1));
    
    a=0.01;
    b=1;
    M(5,e)=(a+b)/2;
    for j=1:50
        LHS=t0(6,e)/t0(5,e);
        RHS=(1+(gamma_a-1)/2*M(6,e)*M(6,e))/(1+(gamma_a-1)/2*M(5,e)*M(5,e))*(M(6,e)*(1+gamma_a*M(5,e)*M(5,e))/M(5,e)/(1+gamma_a*M(6,e)*M(6,e)))^(2);
        if(LHS-RHS)<0
            a=M(5,e);
            M(5,e)=(a+b)/2;
        end
        if(LHS-RHS)>0
            b=M(5,e);
            M(5,e)=(a+b)/2;
        end
        if LHS==RHS
            M(5,e)=(a+b)/2;
            break;
        end
    end
    
    t(5,e)=t0(5,e)/(1+(gamma_a-1)/2*M(5,e)*M(5,e));
    t(6,e)=t0(6,e)/(1+(gamma_a-1)/2*M(6,e)*M(6,e));
    
    A_comb(1,e)=A_crit(1,e)*(1/M(5,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(5,e)*M(5,e)))^((gamma_a+1)/2/(gamma_a-1)));
    v(4,e)=M(4,e)*(gamma_a*R_a*t(4,e))^(1/2);
    v(5,e)=M(5,e)*(gamma_a*R_a*t(5,e))^(1/2);
    rho(5,e)=rho(4,e)*A_comb_in*v(4,e)/v(5,e)/A_comb(1,e);
    p(5,e)=rho(5,e)*R_a*t(5,e);
    p(6,e)=p(5,e)*(1+gamma_a*M(5,e)*M(5,e))/(1+gamma_a*M(6,e)*M(6,e));
    p0(5,e)=p(5,e)*(t0(5,e)/t(5,e))^(gamma_a+1/gamma_a-1);
    p0(6,e)=p(6,e)*(t0(6,e)/t(6,e))^(gamma_a+1/gamma_a-1);
    
    M(7,e)=M(6,e);
    p(7,e)=p(6,e);
    t(7,e)=t(6,e);
    p0(7,e)=p0(6,e);
    t0(7,e)=t0(6,e);
    rho(7,e)=rho(6,e);
    
    v(7,e)=M(7,e)*(gamma_a*R_a*t(7,e))^(1/2);
    
    %Nozzle
    
    A_crit(3,e)=A_comb(1,e)/(1/M(7,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(7,e)*M(7,e)))^((gamma_a+1)/2/(gamma_a-1)));
    A_exit=he*w; 

    guess_Noz3=[1 30];

    M(8,e)=fzero(@Nozzle,guess_Noz3,[],A_exit,A_crit(3,e)); 
    t(8,e)=t(7,e)*(cp + gamma_a*M(7,e)*M(7,e)*R_a/2)/(cp + gamma_a*M(8,e)*M(8,e)*R_a/2);
    p(8,e)=p(7,e)*(t(8,e)/t(7,e))^(gamma_a/(gamma_a-1));
    p0(8,e)=p(8,e)*(1+(gamma_a-1)/2*M(8,e)*M(8,e))^(gamma_a/(gamma_a-1));
    t0(8,e)=t0(7,e);

    v(8,e)=M(8,e)*(gamma_a*R_a*t(8,e))^(1/2);
    v(8,e)=M(8,e)*(gamma_a*R_a*t(8,e))^(1/2);

    net_thrust(1,e)=(mdot(1,e)+mdot_fuel(1,e))*v(8,e) - mdot(1,e)*v(1,e) + p(8,e)*he*w - p(1,e)*hc*w;
    nozzle_thrust(1,e)=(mdot(1,e)+mdot_fuel(1,e))*(v(8,e) - v(7,e)) + p(8,e)*he*w - p(7,e)*h3*w;%Pressure Thrust
    g0=9.812;
    isp(1,e)=(net_thrust(1,e))/g0/mdot_fuel(1,e);
    
    t_eff(1,e)=((mdot(1,e)+mdot_fuel(1,e))*v(8,e)*v(8,e)/2 - mdot(1,e)*v(1,e)*v(1,e)/2)/mdot_fuel(1,e)/lhv;
    p_eff(1,e)=2*v(1,e)/v(8,e)/(1+v(1,e)/v(8,e));
    o_eff(1,e)=t_eff(1,e)*p_eff(1,e);
end


   end
        if LHS==RHS
            M(5,e)=(a+b)/2;
            break;
        end
    end
    
    t(5,e)=t0(5,e)/(1+(gamma_a-1)/2*M(5,e)*M(5,e));
    t(6,e)=t0(6,e)/(1+(gamma_a-1)/2*M(6,e)*M(6,e));
    
    A_comb(1,e)=A_crit(1,e)*(1/M(5,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(5,e)*M(5,e)))^((gamma_a+1)/2/(gamma_a-1)));
    v(4,e)=M(4,e)*(gamma_a*R_a*t(4,e))^(1/2);
    v(5,e)=M(5,e)*(gamma_a*R_a*t(5,e))^(1/2);
    rho(5,e)=rho(4,e)*A_comb_in*v(4,e)/v(5,e)/A_comb(1,e);
    p(5,e)=rho(5,e)*R_a*t(5,e);
    p(6,e)=p(5,e)*(1+gamma_a*M(5,e)*M(5,e))/(1+gamma_a*M(6,e)*M(6,e));
    p0(5,e)=p(5,e)*(t0(5,e)/t(5,e))^(gamma_a+1/gamma_a-1);
    p0(6,e)=p(6,e)*(t0(6,e)/t(6,e))^(gamma_a+1/gamma_a-1);
    
    M(7,e)=M(6,e);
    p(7,e)=p(6,e);
    t(7,e)=t(6,e);
    p0(7,e)=p0(6,e);
    t0(7,e)=t0(6,e);
    rho(7,e)=rho(6,e);
    
    v(7,e)=M(7,e)*(gamma_a*R_a*t(7,e))^(1/2);
    
    %Nozzle
    
    A_crit(3,e)=A_comb(1,e)/(1/M(7,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(7,e)*M(7,e)))^((gamma_a+1)/2/(gamma_a-1)));
    A_exit=he*w; 

    guess_Noz3=[1 30];

    M(8,e)=fzero(@Nozzle,guess_Noz3,[],A_exit,A_crit(3,e)); 
    t(8,e)=t(7,e)*(cp + gamma_a*M(7,e)*M(7,e)*R_a/2)/(cp + gamma_a*M(8,e)*M(8,e)*R_a/2);
    p(8,e)=p(7,e)*(t(8,e)/t(7,e))^(gamma_a/(gamma_a-1));
    p0(8,e)=p(8,e)*(1+(gamma_a-1)/2*M(8,e)*M(8,e))^(gamma_a/(gamma_a-1));
    t0(8,e)=t0(7,e);

    v(8,e)=M(8,e)*(gamma_a*R_a*t(8,e))^(1/2);
    v(8,e)=M(8,e)*(gamma_a*R_a*t(8,e))^(1/2);

    net_thrust(1,e)=(mdot(1,e)+mdot_fuel(1,e))*v(8,e) - mdot(1,e)*v(1,e) + p(8,e)*he*w - p(1,e)*hc*w;
    nozzle_thrust(1,e)=(mdot(1,e)+mdot_fuel(1,e))*(v(8,e) - v(7,e)) + p(8,e)*he*w - p(7,e)*h3*w;%Pressure Thrust
    g0=9.812;
    isp(1,e)=(net_thrust(1,e))/g0/mdot_fuel(1,e);
    
    t_eff(1,e)=((mdot(1,e)+mdot_fuel(1,e))*v(8,e)*v(8,e)/2 - mdot(1,e)*v(1,e)*v(1,e)/2)/mdot_fuel(1,e)/lhv;
    p_eff(1,e)=2*v(1,e)/v(8,e)/(1+v(1,e)/v(8,e));
    o_eff(1,e)=t_eff(1,e)*p_eff(1,e);
end


