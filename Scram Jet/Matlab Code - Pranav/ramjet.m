






%intial conditions
Mdes=9.52;
alt=26; %altitude: 26 km 
R_a=287.1;
R_g=284.09;
p1=2163.2;     %in Pa
t1=219.34; %in K
rho1=p1/R_a/t1;
gamma_a=1.4;
gamma_g=1.33;
cp=gamma_a*R_a/(gamma_a-1);

%body parameters
h1=0.2855;
h2=0.63;
h3=0.0845;
hc=h1+h2+h3; %capture height
he=1.5;      %exit height
w=1;         %width of engine

beta=zeros(3,21);
theta=zeros(3,21);

phi=0.3;

%defining fluid property matrices
M=zeros(8,21);
v=zeros(8,21);
p=zeros(8,21);
t=zeros(8,21);
p0=zeros(8,21);
t0=zeros(8,21);
rho=zeros(8,21);

%defining arrays for different parameters
inlet_drag=zeros(1,21);
inlet_lift=zeros(1,21);
totalpressureloss=zeros(1,21);
A_crit=zeros(3,21);
nozzle_thrust=zeros(1,21);
net_thrust=zeros(1,21);
isp=zeros(1,21);

t_eff=zeros(1,21);
p_eff=zeros(1,21);
o_eff=zeros(1,21);

for i = 1:21
    theta(1,i)=4/180*pi;
    theta(2,i)=8/180*pi;
    theta(3,i)=12/180*pi;   
    
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
        guess=[85*pi/180 asin(1/M(i,e))+pi/180];
    
        beta(i,e)=fzero(@ThetaBetaSolve,guess,[],theta(i,e),M(i,e));
    
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
    
    A_comb=h3*w;
    A_crit(1,e)=A_comb/(1/M(4,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(4,e)*M(4,e)))^((gamma_a+1)/2/(gamma_a-1)));
    
    guess_Noz1=[0.01 1];

    M(5,e)=fzero(@Nozzle,guess_Noz1,[],A_comb,A_crit(1,e)); 
    t(5,e)=t(4,e)*(cp + gamma_a*M(4,e)*M(4,e)*R_a/2)/(cp + gamma_a*M(5,e)*M(5,e)*R_a/2);
    p(5,e)=p(4,e)*(t(5,e)/t(4,e))^(gamma_a/(gamma_a-1));
    p0(5,e)=p(5,e)*(1+(gamma_a-1)/2*M(5,e)*M(5,e))^(gamma_a/(gamma_a-1));
    t0(5,e)=t0(4,e);
    
    t0(6,e)=t0(5,e) + q(1,e)/(gamma_a*R_a/(gamma_a-1));

    t0max=t0(5,e)/(2*(1+gamma_a)*M(5,e)*M(5,e))/(1+(gamma_a-1)/2*M(5,e)*M(5,e))*(1+gamma_a*M(5,e)*M(5,e))^(2);
    
    a=1;
    b=M(5,e);
    M(6,e)=(a+b)/2;
    for j=1:50
        LHS=t0max/t0(6,e);
        RHS=1/(2*(1+gamma_g)*M(6,e)*M(6,e))/(1+(gamma_g-1)/2*M(6,e)*M(6,e))*(1+gamma_g*M(6,e)*M(6,e))^(2);
        if(LHS-RHS)>0
            a=M(6,e);
            M(6,e)=(a+b)/2;
        end
        if(LHS-RHS)<0
            b=M(6,e);
            M(6,e)=(a+b)/2;
        end
        if LHS==RHS
            M(6,e)=(a+b)/2;
            break;
        end
    end
    
    p(6,e)=p(5,e)*(1+gamma_a*M(5,e)*M(5,e))/(1+gamma_a*M(6,e)*M(6,e));
    t(6,e)=t0(6,e)/(1+(gamma_a-1)/2*M(6,e)*M(6,e));
    p0(6,e)=p(6,e)*(1+(gamma_a-1)/2*M(6,e)*M(6,e))^(gamma_a/(gamma_a-1));

    %Diffuser 2
    
    A_comb=h3*w;
    A_crit(2,e)=A_comb/(1/M(6,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(6,e)*M(6,e)))^((gamma_a+1)/2/(gamma_a-1)));
    
    guess_Noz1=[1 30];

    M(7,e)=fzero(@Nozzle,guess_Noz1,[],A_comb,A_crit(2,e)); 
    t(7,e)=t(6,e)*(cp + gamma_a*M(6,e)*M(6,e)*R_a/2)/(cp + gamma_a*M(7,e)*M(7,e)*R_a/2);
    p(7,e)=p(6,e)*(t(7,e)/t(6,e))^(gamma_a/(gamma_a-1));
    p0(7,e)=p(7,e)*(1+(gamma_a-1)/2*M(7,e)*M(7,e))^(gamma_a/(gamma_a-1));
    t0(7,e)=t0(6,e);

    %Nozzle
    
    A_crit(3,e)=A_comb/(1/M(7,e)*(2/(gamma_a+1)*(1+(gamma_a-1)/2*M(7,e)*M(7,e)))^((gamma_a+1)/2/(gamma_a-1)));
    A_exit=he*w; 

    guess_Noz3=[1 30];

    M(8,e)=fzero(@Nozzle,guess_Noz3,[],A_exit,A_crit(3,e)); 
    t(8,e)=t(7,e)*(cp + gamma_a*M(7,e)*M(7,e)*R_a/2)/(cp + gamma_a*M(8,e)*M(8,e)*R_a/2);
    p(8,e)=p(7,e)*(t(8,e)/t(7,e))^(gamma_a/(gamma_a-1));
    p0(8,e)=p(8,e)*(1+(gamma_a-1)/2*M(8,e)*M(8,e))^(gamma_a/(gamma_a-1));
    t0(8,e)=t0(7,e);

    v(8,e)=M(8,e)*(gamma_a*R_a*t(8,e))^(1/2);
    v(7,e)=M(7,e)*(gamma_a*R_a*t(7,e))^(1/2);

    net_thrust(1,e)=(mdot(1,e)+mdot_fuel(1,e))*v(8,e) - mdot(1,e)*v(1,e) + p(8,e)*he*w - p(1,e)*hc*w;
    nozzle_thrust(1,e)=(mdot(1,e)+mdot_fuel(1,e))*(v(8,e) - v(7,e)) + p(8,e)*he*w - p(7,e)*h3*w;%Pressure Thrust
    g0=9.812;
    isp(1,e)=(net_thrust(1,e))/g0/mdot_fuel(1,e);
    
    t_eff(1,e)=((mdot(1,e)+mdot_fuel(1,e))*v(8,e)*v(8,e)/2 - mdot(1,e)*v(1,e)*v(1,e)/2)/mdot_fuel(1,e)/lhv;
    p_eff(1,e)=2*v(1,e)/v(8,e)/(1+v(1,e)/v(8,e));
    o_eff(1,e)=t_eff(1,e)*p_eff(1,e)
end

