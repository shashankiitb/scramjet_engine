% Creating Function getting M2 after fanno flow.
function fannoflowintergate = fannoint(M1,M2) 
    D = 0.4;
    L = 0.5;
    gamma = 1.4;
    f = 2*10^-3;
    fannoflowintergate = -4*f*L/D + 1/gamma/M1^2 + ((gamma+1)/(2*gamma))*log((M1^2)/((1+(gamma-1)/2*M1^2)))
    - (1/gamma/M2^2 + ((gamma+1)/(2*gamma))*log((M2^2)/((1+(gamma-1)/2*M2^2))));
end 