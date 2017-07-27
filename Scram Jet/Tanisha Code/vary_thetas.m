%Program to give outputs at each station of inlet for multiple combinations
%of theta.

%You can give mach number and range of theta in the inputs section of code

%Temperatures and pressures at station 4 are displayed. However, all
%gas properties are stored in the workspace.

%Graph of variation of stagnation pressure is plotted, can be changed to
%temperature.

%The program also gives you possible thetas for a desired output condition.
%If you wish to modify the particular output condition, see line 114.Here I have
%found condition for maximum stagnation pressure which turns out to be the
%minimum theta values.

%careful of syntax case - T TEMPERTAURES and t thetas  
clc;clear;  
%INPUTS
    Min=6;      % Inlet M
    theta1=(7:11)'*pi/180;  %theta range
    theta2=(8:12)'*pi/180;
    p1=2188;         %in Pa
    T1=222.5;         %in K 
    g=1.4; 
    R_a=287.1;  
    rho1=p1/R_a/T1;
    
    %display inputs
    disp(['THETA 1 RANGE :',num2str((theta1')*180/pi)]);
    disp(['THETA 2 RANGE :',num2str((theta2')*180/pi)]);
    disp(['Inlet M  :',num2str(Min)]);
    disp(['P (Pa)   :',num2str(p1)]);
    disp(['T (K)    :',num2str(T1)]);
    disp(['gamma    :',num2str(g)]);
    disp(['Density  :',num2str(rho1)]);
      
  %define variables
    %range of theta
    t1R=max(size(theta1));
    t2R=max(size(theta2));  
    %Numbers refer to station position
    M4=zeros(t1R,t2R);M3=zeros(t1R,t2R);M2=zeros(t1R,t2R);
    beta3=zeros(t1R,t2R);beta2=zeros(t1R,t2R);beta1=zeros(t1R,t2R);
    T2=zeros(t1R,t2R);T02=zeros(t1R,t2R);  
    P2=zeros(t1R,t2R);P02=zeros(t1R,t2R);
    T3=zeros(t1R,t2R);T03=zeros(t1R,t2R);  
    P3=zeros(t1R,t2R);P03=zeros(t1R,t2R);
    T4=zeros(t1R,t2R);T04=zeros(t1R,t2R);  
    P4=zeros(t1R,t2R);P04=zeros(t1R,t2R);
    
    %Variables for loop
    M=zeros(4,1);
    beta=zeros(3,1);
    start=zeros(4,1); %for beta
    v=zeros(4,1);        %Velocity
    P=zeros(4,1);        %Pressure
    T=zeros(4,1);        %Temperature
    P0=zeros(4,1);    %Total Pressure
    T0=zeros(4,1);    %Total Temperature
    rho=zeros(4,1);    %Density
   
    
%defining initial conditions
    M(1,1)=Min;    
    P(1,1)=p1;
    T(1,1)=T1;
    rho(1,1)=rho1;
    v(1,1)=M(1,1)*(g*R_a*T(1,1))^(1/2);
    T0(1,1)=T(1,1)*(1+0.2*M(1,1)*M(1,1));
    P0(1,1)=P(1,1)*(T0(1,1)/T(1,1))^(g/(g-1));

% THETA LOOP 
  for t1=1:t1R
      for t2=1:t2R
          theta=[theta1(t1,1);theta2(t2,1);theta1(t1,1)+theta2(t2,1)]; 
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
          
            x=(1:4)';
            plot(x,P0,'color',[1/t1 0 1/t2]);  
            xlabel('Station');
            ylabel('Stagnation Pressure')
            legend(['theta1 :',num2str(theta1(t1,1)*180/pi),'  theta2 :',num2str(theta2(t2,1)*180/pi)]);
            hold on;
            pause(0.1);
          
          beta1(t1,t2)=beta(1,1)*180/pi;
          beta2(t1,t2)=beta(2,1)*180/pi;
          beta3(t1,t2)=beta(3,1)*180/pi;
          M2(t1,t2)=M(2,1);M3(t1,t2)=M(3,1);M4(t1,t2)=M(4,1);
          T2(t1,t2)=T(2,1);T3(t1,t2)=T(3,1);T4(t1,t2)=T(4,1);
          P2(t1,t2)=P(2,1);P3(t1,t2)=P(3,1);P4(t1,t2)=P(4,1);
          T02(t1,t2)=T0(2,1);T03(t1,t2)=T0(3,1);T04(t1,t2)=T0(4,1);
          P02(t1,t2)=P0(2,1);P03(t1,t2)=P0(3,1);P04(t1,t2)=P0(4,1);
      end
      
      if t1~=t1R
      pause(0.5);
      legend('CHANGING THETA 1')
      hold on;
      end
  end
  
  %Edit this search function as per requirement. eg. T>800K
  [row,column]=find(P04==max(max(P04)));   
  disp('For minimum stagnation pressure loss :');
  disp(['Theta 1 = ',num2str(theta1(row,1)*180/pi)]);
  disp(['Theta 2 = ',num2str(theta2(column,1)*180/pi)]);
  
  disp(' ');
  disp('Columns : THETA 2    Rows : THETA 1');
  disp('Temperatures at combustor inlet');
  A=[theta1*180/pi T4];
  B=[0 (theta2')*180/pi];
  CT=[B;A];
  disp(CT);
  
  disp('Pressures at combustor inlet');
  A=[theta1*180/pi P4];
  B=[0 (theta2')*180/pi];
  CP=[B;A];
  disp(CP);