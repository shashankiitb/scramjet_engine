% For Fanno Flow Analysis.
% Author : Shashank Verma
% Email : shashankoxy@gmail.com and shashankverma@iitb.ac.in



function fanno_flow_2()
    global gamma           % Declearing Global Variable gamma.
    % Inlet Condition.
    M1 = 3.5; %Inlet Mach no.
    L  = 28;  %Length of Pipe in meters.
    D  = 0.4;   %Diameter of Pipe in meeters.
    f  = 2*10^-3; % Coefficent of Friction.
    gamma = 1.4; % Gamma of air.
    choke = 0; % Chocking of flow condition, if 0 unchocked, if 1 chocked.\
    Lchoke = 0;

    x1 = 0; % Starting Point of combustor.
    x2 = L; % End Point of Combustor.

    %   Checking for chocked flow i.e L* >= L.

    Lchoke=fzero(@fannoflowchokelength,0.0,[],M1);
    disp(Lchoke);% Lchoke is length from start of combustor.
    global fm2
    if (Lchoke - L)>0   % i.e chocking doesn't happen
        p1m = ppm(M1);
        fm1 = ffm(M1);
        fm2 = fm1 - 4*f*L/D;
        M2  = fsolve(@ffmi,M1)
    end
    
    
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defining Pressure ratio function i.e P/P* = 1/M*((gamma+1)/(2+(gamma+1)*M^2))^0.5
function y = ppm(M)
    global gamma
    y = 1/M*((gamma+1)/(2+(gamma+1)*M^2))^0.5;
end

% Defining Friction function i.e 4*f*l/D = (1-M^2)/gamma/M^2 + (gamma+1)/2/gamma* log(((gamma+1)*M^2)/2/(2+(gamma-1)*M^2))
function y = ffm(M)
    global gamma
    y = (1-M^2)/gamma/M^2 + (gamma+1)/2/gamma* log(((gamma+1)*M^2)/(2+(gamma-1)*M^2));
end

% Defining inverse friction function i.e to get value of M from this term value 4*f*l/D
function y = ffmi(M)
    global fm2
    global gamma
    y = -fm2 + (1-M^2)/gamma/M^2 + (gamma+1)/2/gamma* log(((gamma+1)*M^2)/(2+(gamma-1)*M^2));
end   
    
