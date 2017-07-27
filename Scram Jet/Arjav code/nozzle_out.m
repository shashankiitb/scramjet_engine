function [nozzle_Thrust,M6,T6,P6]=nozzle_out(M5,P5,T5,ga,Ar,A_in)
aastrin=((((ga+1)/2)^((-ga-1)/(2*(ga-1))))*((1+((ga-1)/2)*M5*M5)^((ga+1)/(2*(ga-1)))))/M5;
aastrout=aastrin*Ar;
syms M;
eqn=[aastrout==((((ga+1)/2)^((-ga-1)/(2*(ga-1))))*((1+((ga-1)/2)*M*M)^((ga+1)/(2*(ga-1)))))/M];
m6=solve(eqn,M);
M6r=double(m6);
x=find(imag(M6r)==0);
M6p(:,1)=M6r(x(:,1));
M6=M6p(M6p>1);
x=(1+((ga-1)/2)*M5*M5);
t05=T5*x;
T6=t05/(1+((ga-1)/2)*M6*M6);
p05=P5*(x^(ga/(ga-1)));
P6=p05/((1+((ga-1)/2)*M6*M6)^(ga/(ga-1)));
MVin=(P5/((287)*T5))*A_in*(M5*((ga*287*T5)^(0.5)));
MVout=(P6/((287)*T6))*A_in*(M6*((ga*287*T6)^(0.5)));
nozzle_Thrust=(MVout-MVin)+(P6*(Ar*A_in)-P5*A_in);
end