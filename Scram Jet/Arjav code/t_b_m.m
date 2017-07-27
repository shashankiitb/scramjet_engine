function betdg=t_b_m(M1,the1,ga)
syms be
relat=[(tan(the1))==(2*cot(be))*(((M1*M1*sin(be)*sin(be))-1)/((M1*M1)*(ga+cos(2*be))+2))];
be=solve(relat,be);
doubbe=real(double(be));
pos=find(doubbe>0);
c(:,1)=doubbe(pos(:,1));
bet=min(c);
betdg=bet*180/pi;
end