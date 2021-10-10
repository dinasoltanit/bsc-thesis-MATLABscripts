clc; clear; close all;
%% IHT of the GF
lex=0.01; lem=0.01; gap=0; epsr=2.7; epsp=1-25*1j; eps0=8.85*1e-12; h=0.17; %d=0.003; 
c1=0; c2=0; c3=0; y=1;
taw = [0.01689, 0.2206];
A = [0.06166, 0.3204];
n=2; kmax=2*n-1; 
N= 10*2;
Z1=zeros(N,N); b1=zeros(N,1);
dx=0.1;
x1= [0:dx:1 1+dx:dx:2]; z1= [-h,0]; %mm  !!! >>>> gap is zero
%x1n=0:dx:10; z1n=0.17:dz:0.23; %mm
%% Functions
%coeffs. #1
A1= @(r,d) (exp(-r.*d)).*((epsp* cosh(r.*d) .* sinh(r.*h) + epsr .* cosh(r.*h) .* sinh(r.*d))./((2*pi.*r).*(epsr* cosh(r.*h) + epsp* sinh(r.*h))));
%p1=A1(1)
fe11= @(r,d,zm) A1(r,d).* exp(-r.*zm + r.*d +1);
%p2=fe11(1)
%f11
Fe11=@(r,d,zm) (fe11(r,d,zm)) - double(Serii(fe11,0.0001,kmax,d,zm)); %new
%p3=Fe11(1)
F1e11=@(r,w,d,zm) (Fe11(r,d,zm)).*(r.^(-0.5)).*(exp(-1i*w.*r)).*(whittakerW(0,0,-2*1i*w.*r)).*exp(-r);
Saay1= @(t,w,d,zm) (exp(-pi*1i/4))*(F1e11(1i*t./w,w,d,zm)) - (exp(pi*1i/4))*(F1e11(-1i*t./w,w,d,zm));
saay_sum1=@(w,d,zm) A.*Saay1(taw,w,d,zm);
saay_wt1=@(w,d,zm) sum(saay_sum1(w,d,zm));
if11=@(w,d,zm) (((1i/w)./(sqrt(2*pi.*w))).* saay_wt1(w,d,zm)) + double(Serii(fe11,0.0001,kmax,d,zm)).*(F1to2(kmax, w));

%coeffs. #2
B =@(r,d) (exp(-r.*d)).*((epsp * sinh(r.*h) - epsr * cosh(r.*h))./((4*pi.*r).*(epsr* cosh(r.*h) + epsp* sinh(r.*h))));
C =@(r,d) (exp(-r.*d))./(4*pi.*r);
V12=@(r,d,zm) B(r,d) .* exp(-r.*zm) + C(r,d).* exp(r.*zm); % for z_prime > z > 0
fe12=@(r,d,zm) (exp(r)).*V12(r,d,zm);
%f12
Fe12=@(r,d,zm) (fe12(r,d,zm)) - double(Serii(fe12,0.0001,kmax,d,zm)); %new
F1e12=@(r,w,d,zm) (Fe12(r,d,zm)).*(r.^(-0.5)).*(exp(-1i*w.*r)).*(whittakerW(0,0,-2*1i*w.*r)).*exp(-r);
Saay2=@(t,w,d,zm) (exp(-pi*1i/4))*((F1e12(1i*t./w,w,d,zm))) - (exp(pi*1i/4))*((F1e12(-1i*t./w,w,d,zm)));
saay_sum2=@(w,d,zm) A.*Saay2(taw,w,d,zm);
saay_wt2=@(w,d,zm) sum(saay_sum2(w,d,zm));
if12=@(w,d,zm) (((1i./w)./(sqrt(2*pi*w))).* saay_wt2(w,d,zm)) + double(Serii(fe12,0.0001,kmax,d,zm)).*(F1to2(kmax, w));

%coeffs. #3
D =@(r,d) (-epsp * exp(-r.*d))./((4*pi.*r).*(epsr* cosh(r.*h) + epsp* sinh(r.*h)));
V13=@(r,d,zm) D(r,d) .* (exp(-r.*(zm+h)) - exp(r.*(zm+h))); % for 0 > z > -h
fe13=@(r,d,zm) (exp(r)).*V13(r,d,zm);
%f13
Fe13=@(r,d,zm) (fe13(r,d,zm)) - double(Serii(fe13,0.0001,kmax,d,zm)); %new %inverted sym to double
F1e13=@(r,w,d,zm) (Fe13(r,d,zm)).*(r.^-0.5).*(exp(-1i*w.*r)).*(whittakerW(0,0,-2*1i*w.*r)).*exp(-r);
Saay3=@(t,w,d,zm) (exp(-pi*1i/4))*((F1e13(1i*t./w,w,d,zm))) - (exp(pi*1i/4))*((F1e13(-1i*t./w,w,d,zm)));
saay_sum3=@(w,d,zm) A.*Saay3(taw,w,d,zm);
saay_wt3=@(w,d,zm) sum(saay_sum3(w,d,zm));
if13 =@(w,d,zm) (((1i./w)./(sqrt(2*pi*w))).* saay_wt3(w,d,zm)) + double(Serii(fe13,0.0001,kmax,d,zm)).*(F1to2(kmax, w));

%% -- first eq.
for m1=1:N
    xm1= 0.05+0.1*mod(n-1,10);
    zm1= 0.005+0.01*mod(n-1,6);
   for n1=1:N
       xn1= 0.05+0.1*mod(n-1,10);
       zn1= 0.005+0.01*mod(n-1,6);
      if (zn1 <= zm1)
            fprintf("in 1st domain. \n")
            fprintf("m= %d \n", m1)
            fprintf("n= %d \n", n1)
            c1=c1+1;
            % 1
            tic
            %saay_wt1=0;
            d=zn1;
            zm = zm1;
            % Z matrix
            rho = sqrt(xm1^2+1);
            Z1(m1,n1)= dx*dz*(double(if11(rho,d,zm)));
            fprintf("Zmn= %d \n", Z1(m1,n1))
            fprintf("counter1= %d \n", c1)
            toc
      elseif ( (0 <= zm1)&& (zm1 < zn1) )
            % 2
            fprintf("in 2nd domain. \n")
            fprintf("m= %d \n", m1)
            fprintf("n= %d \n", n1)
            c2=c2+1;
            %saay_wt2=0;
            d=zn1;
            zm = zm1;
            rho = sqrt(xm1^2+1);
            Z1(m1,n1)= dx*dz*(double(if12(rho,d,zm)));
            fprintf("Zmn= %d \n", Z1(m1,n1))
            fprintf("counter1= %d \n", c2)
      else
          fprintf("An error has occured. \n")
      end
   end
end
b1(:,1)= -Vt*eps0;
a1= b1/Z1;
%% -- Second eq.
Z2=zeros(N,N); b2=zeros(N,1);
for m2=1:N
    xm2= 10+0.05+0.1*mod(n-1,100);
    zm2= -h+0.005+0.01*mod(n-1,6);
   for n2=1:N
       xn2= 10+0.05+0.1*mod(n-1,100);
       zn2= -h+0.005+0.01*mod(n-1,6);
      if ( -h < zm2 < 0 )
            %3
            fprintf("in 3rd domain. \n")
            fprintf("m= %d \n", m2)
            fprintf("n= %d \n", n2)
            c3=c3+1;
            %saay_wt3=0;
            d=zn2;
            zm = zm2;
            rho = sqrt(xm2^2+1);
            Z2(m2,n2)= dx*dz*(double((if13(rho,d,zm) )));
            fprintf("Zmn= %d \n", Z2(m2,n2))
            fprintf("counter1= %d \n", c3)
      else
          fprintf("An error has occured. \n")
      end
   end
end
a2= b2/Z2;





