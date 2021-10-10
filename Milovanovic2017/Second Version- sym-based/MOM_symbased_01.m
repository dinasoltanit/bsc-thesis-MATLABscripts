clc; clear; close all;
%% IHT of the GF
syms r w t Vt
lex=0.01; lem=0.01; gap=0; epsr=2.7; epsp=1-25*1j; eps0=8.85*1e-12; h=0.17; %d=0.003; 
c1=0; c2=0; c3=0; y=1;
taw = [0.01689, 0.2206];
A = [0.06166, 0.3204];
n=2; kmax=2*n-1; 
N= 100*6;
Z1=zeros(N,N); b1=zeros(N,1);
dx=0.1; dz=0.01; %mm
x1= [0:dx:10 10+dx:dx:20]; z1= [ -h:dz:(-h+0.06)  0:dz:0.06 ]; %mm  !!! >>>> gap is zero
%x1n=0:dx:10; z1n=0.17:dz:0.23; %mm
%% -- first eq.
for m1=1:N
    xm1= 0.05+0.1*mod(n-1,100);
    zm1= 0.005+0.01*mod(n-1,6);
   for n1=1:N
       xn1= 0.05+0.1*mod(n-1,100);
       zn1= 0.005+0.01*mod(n-1,6);
      if (zn1 <= zm1)
            fprintf("in 1st domain. \n")
            fprintf("m= %d \n", m1)
            fprintf("n= %d \n", n1)
            c1=c1+1;
            % 1
            tic
            saay_wt1=0;
            d=zn1;
            % finding coeffs.
            A1(r)= (exp(-r*d))*((epsp* cosh(r*d) * sinh(r*h) + epsr * cosh(r*h) * sinh(r*d))/((2*pi*r)*(epsr* cosh(r*h) + epsp* sinh(r*h))));
            % finding the green's function in freq domain
            V11(r)= A1(r)* exp(-r*(zm1-d));% for z > z_prime
            fe11(r)= (exp(r))*V11(r);
            %f11
            Fe11(r)= (fe11(r)) - double(Serii(@(r)fe11(r),0.0001,kmax)); %new
            F1e11(r,w)= (Fe11(r))*(r^-0.5)*(exp(-1i*w*r))*(whittakerW(0,0,-2*1i*w*r))*exp(-r);
            Saay1(t,w) = (exp(-pi*1i/4))*(subs(F1e11(r,w),r,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e11(r,w),r,-1i*t/w));
            for ii=1:n
              saay_wt1 = saay_wt1 + A(ii)*subs(Saay1(t,w),t,taw(ii));
            end
            if11(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt1) + Serii(@(r)fe11(r),0.0001,kmax)*(F1to2(kmax, w));
            % Z matrix
            rho = sqrt(xm1^2+1);
            %answer= double(subs(if11(w), w, rho ));
            Z1(m1,n1)= dx*dz*(double(subs(if11(w), w, rho )));
            fprintf("Zmn= %d \n", Z1(m1,n1))
            fprintf("counter1= %d \n", c1)
            toc
      elseif ( (0 <= zm1)&& (zm1 < zn1) )
            % 2
            fprintf("in 2nd domain. \n")
            fprintf("m= %d \n", m1)
            fprintf("n= %d \n", n1)
            c2=c2+1;
            saay_wt2=0;
            d=zn1;
            B(r) = (exp(-r*d))*((epsp * sinh(r*h) - epsr * cosh(r*h))/((4*pi*r)*(epsr* cosh(r*h) + epsp* sinh(r*h))));
            C(r) = (exp(-r*d))/(4*pi*r);
            V12(r)= B(r) * exp(-r*zm1) + C(r)* exp(r*zm1); % for z_prime > z > 0
            fe12(r)= (exp(r))*V12(r);
            %f12
            Fe12(r)= (fe12(r)) - double(Serii(@(r)fe12(r),0.0001,kmax)); %new
            F1e12(r,w)= (Fe12(r))*(r^-0.5)*(exp(-1i*w*r))*(whittakerW(0,0,-2*1i*w*r))*exp(-r);
            Saay2(t,w) = (exp(-pi*1i/4))*(subs(F1e12(r,w),r,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e12(r,w),r,-1i*t/w));
            for ii=1:n
              saay_wt2 = saay_wt2 + A(ii)*( subs(Saay2(t,w),t,taw(ii)) );
            end
            if12(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt2) + Serii(@(r)fe12(r),0.0001,kmax)*(F1to2(kmax, w));
            rho = sqrt(xm1^2+1);
            Z1(m1,n1)= dx*dz*(double(subs(if12(w), w, rho )));
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
            saay_wt3=0;
            d=zn2;
            D(r) = (-epsp * exp(-r*d))/((4*pi*r)*(epsr* cosh(r*h) + epsp* sinh(r*h)));
            V13(r)= D(r) * (exp(-r*(zm2+h)) - exp(r*(zm2+h))); % for 0 > z > -h
            fe13(r)= (exp(r))*V13(r);
            %f13
            Fe13(r)= (fe13(r)) - double(Serii(@(r)fe13(r),0.0001,kmax)); %new %inverted sym to double
            F1e13(r,w)= (Fe13(r))*(r^-0.5)*(exp(-1i*w*r))*(whittakerW(0,0,-2*1i*w*r))*exp(-r);
            Saay3(t,w)= (exp(-pi*1i/4))*(subs(F1e13(r,w),r,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e13(r,w),r,-1i*t/w));
            for ii=1:n
              saay_wt3 = saay_wt3 + A(ii)*subs(Saay3(t,w),t,taw(ii));
            end
            if13(w) =(((1i/w)/(sqrt(2*pi*w)))* saay_wt3) + Serii(@(r)fe13,0.0001,kmax)*(F1to2(kmax, w));
            rho = sqrt(xm2^2+1);
            Z2(m2,n2)= dx*dz*(double(subs(if13(w), w, rho )));
            fprintf("Zmn= %d \n", Z2(m2,n2))
            fprintf("counter1= %d \n", c3)
      else
          fprintf("An error has occured. \n")
      end
   end
end
a2= b2/Z2;
% It is done. the only thing remained is to find rho. however, the solution
% is quit rough.
% I believe, the results needs some refinements. 1st, it has to be improved
% in terms of performance. 2nd, the solution can be changed to become based
% on Galerkin's method rather than the point-matching method. page 62/289
% from Gibson





