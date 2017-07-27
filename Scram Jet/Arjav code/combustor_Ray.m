function [P_out,T_out,M_out]=combustor_Ray(M2,T2,P2,ga,cp,q)
q_act=q*1000000;
P2_act=P2*1000;
Tstag_in=T2*(1+((ga-1)/2)*M2*M2);
Tstag_out=(cp*Tstag_in+q_act)/cp;
ttstr4=2*(1+ga)*(((1+((ga-1)/2)*M2*M2)*M2*M2)/((1+ga*M2*M2)*(1+ga*M2*M2)));
ttstr5=((Tstag_out/Tstag_in)*ttstr4);
if ttstr5>1
    disp('Flow is choked within the combustor')
elseif ttstr5==1
    disp('Flow is choked at the end of nozzle so M_exit==1')
    M_out=1;
else
    syms M
    eqn=[ttstr5==2*(1+ga)*(((1+((ga-1)/2)*M*M)*M*M)/((1+ga*M*M)*(1+ga*M*M)))];
    m_out=solve(eqn,M);
    M_outr=double(m_out);
    x=find(M_outr>1);
    M_out=M_outr(x(1));
    T_out=Tstag_out*(1/(1+((ga-1)/2)*M_out*M_out));
    p5p5str=((1+ga)/(1+ga*M_out*M_out));
    p4p4str=((1+ga)/(1+ga*M2*M2));
    p5p4=p5p5str/p4p4str;
    P_out=p5p4*P2;
end
