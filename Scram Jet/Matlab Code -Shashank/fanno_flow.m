% For Fanno Flow Analysis.
% Author : Shashank Verma
% Email : shashankoxy@gmail.com and shashankverma@iitb.ac.in

clc;

% Inlet Condition.
M1 = 3.5; %Inlet Mach no.
L  = 0.5;  %Length of Pipe in meters.
D  = 0.4;   %Diameter of Pipe in meeters.
f  = 2*10^-3; % Coefficent of Friction.
gamma_a = 1.4; % Gamma of air.
choke = 0; % Chocking of flow condition, if 0 unchocked, if 1 chocked.\
Lchoke = 0;

x1 = 0; % Starting Point of combustor.
x2 = L; % End Point of Combustor.

% Checking for chocked flow i.e L* >= L.

Lchoke=fzero(@fannoflowchokelength,0.0,[],M1); % Lchoke is length from start of combustor.

if (Lchoke - L)>0
    guess = [];
    M2 = fzero(@fannoflowintergate,guess,[],M1);
end    
    

    
    