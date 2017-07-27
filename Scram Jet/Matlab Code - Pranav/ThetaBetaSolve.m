function ThetaBetaSolve=thetabetasolve(beta,theta,M)
    gamma=1.4;
    ThetaBetaSolve=cot(theta)-tan(beta)*((gamma+1)*M*M/2/(M*M*sin(beta)*sin(beta)-1)-1);
end 


 