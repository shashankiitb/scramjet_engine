% Below code 


function Nozzle=nozzle(M,A_exit,A_crit)
    gamma=1.4;
    Nozzle=A_exit/A_crit-(1/M*(2/(gamma+1)*(1+(gamma-1)/2*M*M))^((gamma+1)/2/(gamma-1)));
end

