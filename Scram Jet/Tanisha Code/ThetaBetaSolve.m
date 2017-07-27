function ThetaBetaSolve=thetabetasolve(beta,theta,M)
    gamma=1.4;
    b=beta*pi/180;
    tta=theta*pi/180;
    mn=M*sin(b);
    ThetaBetaSolve=tan(b)/tan(b-tta)-(gamma+1)*mn*mn/(2+(gamma-1)*mn*mn);
end 

