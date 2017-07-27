array=input('Provide inlet conditions in this order"[M1 P1 T1 gamma]": ');
format long
M1=array(1,1);
P1=array(1,2);
T1=array(1,3);
ga=array(1,4);
po1 = P1*(1+(ga-1)/2*M1^2)^(ga/ga-1);
i=1;
while 1>0
    the1(i)=input('Please provide the theta here: ');
    betdg=t_b_m(M1,the1(i)*pi/180,ga);
    bet=num2str(betdg);
    x=['The beta here is: ',bet];
    disp(x)
    [M2,P2,T2,Pstag2]=next_m_and_prop(betdg*pi/180,the1(i)*pi/180,M1,ga,P1,T1,po1);
    x=[M2,P2,T2,Pstag2];
    disp('The after-shock values in the order [M,P,T,Pstagnation] are: ')
    disp(x)
    x=input('Input 1 if there are any more inlet angles, 0 if not');
    if x==0
        break
    else
        x=0;
    end
    M1=M2;
    T1=T2;
    P1=P2;
    po1=Pstag2;
    i=i+1;
end
%cowl-
M1=M2;
T1=T2;
P1=P2;
po1=Pstag2;
    
cowltheta=sum(the1);
betdg=t_b_m(M1,cowltheta*pi/180,ga);
bet=num2str(betdg);
x=['The beta at cowl is: ',bet];
disp(x)
[M2,P2,T2,Pstag2]=next_m_and_prop(betdg*pi/180,cowltheta*pi/180,M1,ga,P1,T1,po1);
x=[M2,P2,T2,Pstag2];
disp('The combustor inlet values in the order [M,P,T,Pstagnation] are: ')
disp(x)
q=input('Please input the value of "q" : ');
cp=input('Please input the value of "Cp" : ');
[P5,T5,M5]=combustor_Ray(M2,T2,P2,ga,cp,q);
x=[P5 T5 M5];
disp('The combustor outlet values of [P M T] respectively are: ')
disp(x)
A_in=input('Please input nozzle inlet area: ');
Ar=input('Please input nozzle area ratio (A_out/A_in): ');
[nozzle_Thrust,M6,T6,P6]=nozzle_out(M5,P5,T5,ga,Ar,A_in);
disp('The Values of [nozzle_Thrust M T P] at nozzle exit respectively are: ')
disp([nozzle_Thrust M6 T6 P6])