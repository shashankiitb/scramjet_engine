function [M2,P2,T2,Pstag2]=next_m_and_prop(be,the,M1,ga,p1,t1,po1)
m1=M1*sin(be);
m2=(((m1*m1)+(2/(ga-1)))/(((2*ga*m1*m1)/(ga-1))-1))^(0.5);
M2=m2/(sin(be-the));
p2p1=((1+(ga*m1*m1))/(1+(ga*m2*m2)));
P2=p2p1*p1;
t2t1=((1+(((ga-1)/2)*m1*m1))/(1+(((ga-1)/2)*m2*m2)));
T2=t2t1*t1;
po2po1=(((((ga+1)*m1*m1/2)/(1+((ga-1)*m1*m1)/2))^(ga/(ga-1)))/((p2p1)^(1/(ga-1))));
Pstag2=po2po1*po1;
end