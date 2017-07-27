clc;clear;

gamma = 1.4;
Cf = 0.002;
Dm = 2.75;
xmax = 22.75;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
book_fig_x = [ 0.00369251	0.0585716
0.365559	0.0724021
0.653575	0.0858901
1.00436	0.1122
1.5767	0.132206
2.06042	0.144535
2.51829	0.158339
3.03894	0.168822
3.55589	0.180407
3.95468	0.190193
4.15408	0.20334
4.6858	0.216388
5.18429	0.222842
5.68278	0.2359
6.18127	0.245656
6.67976	0.255412
7.11178	0.261887
7.67674	0.271623
8.17523	0.274775
8.57033	0.28126
8.80665	0.281189
9.17221	0.290985
9.63746	0.0
10.1027	0.303913
10.6012	0.307065
11.0997	0.320123];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pressure_curve_1 = [ 2.28653	0.0542659
9.89296	0.0525794
10.492	0.0525794
10.9684	0.060119
11.4652	0.0606151
11.9639	0.0844246
12.5054	0.115774
13.019	0.137202
13.5045	0.149702
13.9565	0.15744
14.4907	0.170635
14.977	0.176587
15.4967	0.188095
15.9917	0.202381
16.4812	0.221825
16.9966	0.228671
17.5166	0.238393
18.0311	0.251984
18.5756	0.260714
19.052	0.26756
19.5526	0.27629
20.064	0.276587
20.564	0.289782
21.079	0.300496
21.6137	0.30873
22.1054	0.310813
22.606	0.318849
22.8439	0.325496 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pressure_curve_2 = [11.0182	0.0517857
11.588	0.0517857
12.0456	0.0539683
12.479	0.0549603
13.0208	0.0456349
13.5073	0.049504
14.0337	0.0472222
14.5258	0.0459325
15.0044	0.0743056
15.5707	0.102579
16.0159	0.125694
16.4909	0.144246
17.05	0.152679
17.5161	0.164087
18.0797	0.175099
18.5509	0.185119
19.0758	0.194048
19.527	0.20873
20.0567	0.217857
20.5667	0.229365
21.0282	0.239187
21.5874	0.24623
22.034	0.25873
22.5687	0.266865
22.8409	0.271329 ];
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pressure_curve_3 = [14.0523	0.0542659
14.5686	0.0535714
15.0074	0.05
15.465	0.0521825
16.03	0.0524802
16.5411	0.0544643
17.0628	0.0892857
17.5417	0.115179
18.0655	0.133036
18.5895	0.149107
19.0269	0.157143
19.5512	0.170536
20.1199	0.180357
20.5563	0.195536
21.0812	0.204464
21.562	0.216071
22.0423	0.23125
22.5675	0.2375
22.7864	0.240179 ];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pressure_curve_4 = [17.0675	0.0517857
17.5058	0.0517857
18.031	0.0589286
18.0712	0.0875
18.5072	0.10625
19.5556	0.135714
20.0798	0.15
20.5606	0.161607
21.1291	0.172321
21.5224	0.182143
22.091	0.191964
22.5714	0.20625
22.8338	0.211607 ];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pressure_curve_5 = [18.0313	0.05625
18.5144	0.0491071
19.0836	0.0544643
19.5656	0.05625
20.0456	0.0732143
20.5662	0.116964
21.0897	0.136607
21.6139	0.151786
22.0511	0.160714
22.6191	0.175893
22.7938	0.180357 ];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pressure_curve_6 = [18.0215	0.0570437
18.4703	0.0508929
18.9959	0.0544643
19.5659	0.0535714
20.0467	0.0642857
20.6117	0.103571
21.0479	0.120536
21.5722	0.134821
22.0091	0.146429
22.5342	0.153571
22.8401	0.160714 ];

%syms M(x) R(x) P(x)

% Defining Diff Equation   
% x(1) = M(x), x(2) = R(x), x(3) = P(x).
%ode1 = diff(x(1)) == -(x(1)/2)*(1+((gamma-1)*(x(1)*x(1)))/2)*(93*Cf/(Dm*x(2)));
%ode2 = diff(x(2)) == Cf/(2*Dm)*(93+x(1)*x(1)*(93*(gamma-1)-89*gamma*x(2)));
%ode3 = diff(x(3)) == 89/Dm*Cf*gamma*(x(1)*x(1))*x(3)/2;
%odes = [-(x(1)/2)*(1+((gamma-1)*(x(1)*x(1)))/2)*(93*Cf/(Dm*x(2))); Cf/(2*Dm)*(93+x(1)*x(1)*(93*(gamma-1)-89*gamma*x(2))); 89/Dm*Cf*gamma*(x(1)*x(1))*x(3)/2]

f = @(t,x)[-(x(1)/2)*(1+((gamma-1)*(x(1)*x(1)))/2)*(93*Cf/(Dm*x(2))); Cf/(2*Dm)*(93+x(1)*x(1)*(93*(gamma-1)-89*gamma*x(2))); 89/Dm*Cf*gamma*(x(1)*x(1))*x(3)/2];

[t1,xa1] = ode45(f,[0,xmax],[2.60,1.0,2.61]);
[t2,xa2] = ode45(f,[0,xmax],[2.60,1.0,3.103]);
[t3,xa3] = ode45(f,[0,xmax],[2.60,1.0,3.509]);
[t4,xa4] = ode45(f,[0,xmax],[2.60,1.0,4.147]);
[t5,xa5] = ode45(f,[0,xmax],[2.60,1.0,4.727]);
[t6,xa6] = ode45(f,[0,xmax],[2.60,1.0,5.2896]);
a1 = (xa1(:,3))/45.0;
a2 = (xa2(:,3))/53.5;
a3 = (xa3(:,3))/60.5;
a4 = (xa4(:,3))/71.5;
a5 = (xa5(:,3))/81.5;
a6 = (xa6(:,3))/91.2;

for i = 1:45
    if xa1(i,3)<= 14.6959
        na1(i,1)=xa1(i,3);
    end    
end

for i = 1:45
    if xa2(i,3)<= 14.6959
        na2(i,1)=xa2(i,3);
    end    
end

for i = 1:45
    if xa3(i,3)<= 14.6959
        na3(i,1)=xa3(i,3);
    end    
end

for i = 1:45
    if xa4(i,3)<= 14.6959
        na4(i,1)=xa4(i,3);
    end    
end

for i = 1:45
    if xa5(i,3)<= 14.6959
        na5(i,1)=xa5(i,3);
    end    
end

for i = 1:45
    if xa6(i,3)<= 14.6959
        na6(i,1)=xa6(i,3);
    end    
end
plot(t1(1:29)+10.492,na1/45.0);
hold on;
plot(t1(1:23)+14.5254,na2/53.5);
hold on;
plot(t1(1:19)+16.5411,na3/60.5);
hold on;
plot(t1(1:16)+18.031,na4/71.5);
hold on;
plot(t1(1:14)+19.5656,na5/81.5);
hold on;
plot(t1(1:12)+19.5859,na6/91.2);
hold on;
% plot(pressure_curve_1(:,1)-11.9639,pressure_curve_1(:,2),'red');
% hold on;
% plot(pressure_curve_2(:,1)-14.5258,pressure_curve_2(:,2),'blue');
% hold on;
% plot(pressure_curve_3(:,1)-17.0628,pressure_curve_3(:,2),'black');
% hold on;
% plot(pressure_curve_4(:,1)-18.0712,pressure_curve_4(:,2),'green');
% hold on;
% plot(pressure_curve_5(:,1)-20.0456,pressure_curve_5(:,2),'yellow');
% hold on;
% plot(pressure_curve_6(:,1)-20.0467,pressure_curve_6(:,2),'cyan');
% hold off;
plot(pressure_curve_1(:,1),pressure_curve_1(:,2),'.');
hold on;
plot(pressure_curve_2(:,1),pressure_curve_2(:,2),'-.');
hold on;
plot(pressure_curve_3(:,1),pressure_curve_3(:,2),'--');
hold on;
plot(pressure_curve_4(:,1),pressure_curve_4(:,2),'*');
hold on;
plot(pressure_curve_5(:,1),pressure_curve_5(:,2),'+');
hold on;
plot(pressure_curve_6(:,1),pressure_curve_6(:,2),'*-');
hold off;
grid on
grid minor
xlabel('x(m)'), ylabel('Pw');
lgdx = legend(' 45.0 ',' 53.5 ','60.5','71.5','81.5','91.2');%,' 45.0 ',' 53.5 ','60.5','71.5','81.5','91.2');