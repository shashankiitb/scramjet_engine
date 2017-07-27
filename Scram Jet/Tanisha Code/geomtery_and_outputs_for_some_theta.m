clc;clear;      %careful of syntax case - T TEMPERTAURES and t thetas 
   
%INPUTS
    Min=input( 'Enter Mach number  : ');
    tta1=input('Enter theta1 in deg: ')*pi/180;
    tta2=input('Enter theta2 in deg: ')*pi/180;
    p1=2188;         %in Pa
    T1=222.5;         %in K
    g=1.4; 
    R_a=287.1;  
    rho1=p1/R_a/T1;
    
    theta=[tta1;tta2;(tta1+tta2)];
    M=zeros(4,1);
    beta=zeros(3,1);
    start=zeros(4,1);    %for beta
    v=zeros(4,1);        %Velocity
    P=zeros(4,1);        %Pressure
    T=zeros(4,1);        %Temperature
    P0=zeros(4,1);    %Total Pressure
    T0=zeros(4,1);    %Total Temperature
    rho=zeros(4,1);   %Density

    %defining initial conditions
    M(1,1)=Min;
    P(1,1)=p1;
    T(1,1)=T1;
    rho(1,1)=rho1;
    v(1,1)=M(1,1)*(g*R_a*T(1,1))^(1/2);
    T0(1,1)=T(1,1)*(1+0.2*M(1,1)*M(1,1));
    P0(1,1)=P(1,1)*(T0(1,1)/T(1,1))^(g/(g-1));
    
    %Analysing
    for k=1:3       
            start(k,1)=max((asin(1/M(k,1))+pi/18000),(theta(k,1)+pi/18000))*180/pi;
            beta(k,1)=min(fzero(@ThetaBetaSolve,[start(k,1) 45+theta(k,1)/2*180/pi],[],theta(k,1)*180/pi,M(k,1)))*pi/180; 
            
            M(k+1,1) = (((g-1)*M(k,1)*M(k,1)*sin(beta(k,1))*sin(beta(k,1))+2)/(2*g*M(k,1)*M(k,1)*sin(beta(k,1))*sin(beta(k,1))-(g-1))/(sin(beta(k,1)-theta(k,1))*sin(beta(k,1)-theta(k,1))))^(1/2);
            P(k+1,1) = P(k,1)*(2*g*M(k,1)*M(k,1)*sin(beta(k,1))*sin(beta(k,1))-(g-1))/(g+1);
            rho(k+1,1) = rho(k,1)*((g+1)*M(k,1)*M(k,1)*sin(beta(k,1))*sin(beta(k,1)))/((g-1)*M(k,1)*M(k,1)*sin(beta(k,1))*sin(beta(k,1))+2);
            T(k+1,1) = P(k+1,1)/P(k,1)*rho(k,1)/rho(k+1,1)*T(k,1);
            v(k,1) = M(k,1)*(g*R_a*T(k,1))^(1/2);
            T0(k+1,1) = T(k+1,1)*(1+(g-1)/2*M(k+1,1)*M(k+1,1));
            P0(k+1,1) = P(k+1,1)*(T0(k+1,1)/T(k+1,1))^(g/(g-1));                   
    end
      
    % Geomtery point=[y;x]
    hc=input('Enter capture height : ');  
    P1=[hc;0] ;
    xc=hc*cot(beta(1,1));
    Pc=[0;xc]; %cowl
    % Equations y-mx=c
    % y-tan(-theta1)x2=hc  for sfc1
       ms1=tan(-tta1);
       cs1=hc;
    % y-tan(180-(theta1+beta2))x=tan(180-(theta1 +beta2))*(-xc); for shk2
       msk2=tan(pi-(theta(1,1)+beta(2,1)));
       csk2=msk2*(-xc);
    % find intersection of sfc1 and shk2 = P2
       A=[1,-ms1;1,-msk2];
       B=[cs1;csk2];
       P2=A^(-1)*B;
    % SFC2
       ms2=tan(-(tta1+tta2));
       cs2=P2(1,1)-ms2*P2(2,1);
    % SHK3
       msk3=tan(beta(3,1)-(tta1+tta2));
       csk3=msk3*(-xc);
    % intersection of surface 2 and shock3
       A=[1 -ms2;1 -msk3];
       B=[cs2;csk3];
       P3=A^(-1)*B;
              
    % displays
       stn=(1:4)';
       disp(['Station','             M           ','P (Pa)       ','T(K)          ','P0(Pa)        ','T0 (K)']);
       disp([stn,M,P,T,P0,T0]);
       
       disp('Betas of shock 1,2,3');
       disp(beta*180/pi);
          
       disp('POINTS : First row is the y coordinate ');
       disp(['Leading edge     ','First bend  ','Second Bend      ','Cowl']);
       points_XissecondROW=[P1 P2 P3 Pc];
       disp(points_XissecondROW);
         
       plot(points_XissecondROW(2,:),points_XissecondROW(1,:));
       legend(['Mach  ',num2str(Min),'  theta1= ',num2str(theta(1)*180/pi),'  theta2 = ',num2str(theta(2)*180/pi)]);
       disp('Please note, in plot, first two lines are scramjet surfaces and third line is the last shock, thus giving cowl tip position')