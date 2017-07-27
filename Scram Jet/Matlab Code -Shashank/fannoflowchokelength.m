% Creating Function for getting chocked length.
function fannoflowchokelength = fannochoke(Lchoke,M1) 
    D = 0.4;
    gamma = 1.4;
    f = 2*10^-3;
    fannoflowchokelength=-4*f*Lchoke/D + (1-M1^2)/gamma/M1^2 + ((gamma+1)/(2*gamma))*log((M1^2)/((2/(gamma+1))*(1+(gamma-1)/2*M1^2)));
end 
