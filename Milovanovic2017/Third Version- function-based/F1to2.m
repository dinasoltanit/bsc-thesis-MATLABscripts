function [F1to2]= F1to2(m,w)
S1=1;
S2=1;
jj=1;
ii=1;
% temp1=0;
% temp2=0;
err1=10;
err2=10;
while err1 < 0.001
   temp1=S1;
   S1=S1+ ((((-w^2)^(-jj))*(zarb((m + 1)/2, jj))^2)/(factorial(jj) * zarb(0.5, jj)));
   jj=jj+1;
   err1=abs(S1-temp1);  
end
%sum[((((-w^2)^(-j))*(zarb[(m + 1)/2, j])^2)/(factorial(j) * zarb[0.5, j])), {j, 0, inf}] +
while err2 < 0.001
   temp2=S2;
   S2=S2+ ((((-w^2)^(-ii))*(zarb((m + 2)/2, ii))^2)/(factorial(ii) * zarb(1.5, ii)));
   ii=ii+1;
   err2=abs(S2-temp2);  
end
%sum[((((-w^2)^(-j))*(zarb[(m + 2)/2, j])^2)/(factorial(j) * zarb[1.5, j])), {j, 0, inf}];
F1to2 = ((gamma (m + 1)/gamma (1)))*(((gamma(1/2))*(gamma(1))*(w^(-(m + 1))))/((gamma((m + 2)/2))*(gamma((1 - m)/2))))* S1+(((gamma(-1/2))*(gamma(1))*(w^(-(m + 2))))/((gamma((m + 1)/2))*(gamma(-0.5))))*S2;
end