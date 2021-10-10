clc; clear; close all;
%% IHT of the GF
syms r w t x
fs=1e-5; x= fs:fs:0.02;
y=1;
% z had actual values in different sections.
% d had either an actual value
epsr=2.7; epsp=1-25*1j; d=0.002; h=0.001;
z1= -h:1e-4:0;
z2= 1e-4:1e-4:d;
z3= (d+1e-4):1e-4:0.1;

taw = [0.01689, 0.2206, 0.9232, 2.436, 5.125, 9.745];
A = [0.06166, 0.3204, 0.4197, 0.1773, 0.02055, 0.0003617];
n=6;
kmax=2*n-1;
saay_wt1=0;
saay_wt2=0;
saay_wt3=0;

B = (exp(-r*d))*((epsp * sinh(r*h) - epsr * cosh(r*h))/((4*pi*r)*(epsr* cosh(r*h) + epsp* sinh(r*h))));
A1= (exp(-r*d))*((epsp* cosh(r*d) * sinh(r*h) + epsr * cosh(r*h) * sinh(r*d))/((2*pi*r)*(epsr* cosh(r*h) + epsp* sinh(r*h))));
C = (exp(-r*d))/(4*pi*r);
D = (-epsp * exp(-r*d))/((4*pi*r)*(epsr* cosh(r*h) + epsp* sinh(r*h)));

V11(r)= A1* exp(-r*(z1-d)) ;
V12(r)= B * exp(-r*z2) + C* exp(r*z2);
V13(r)= D * (exp(-r*(z3+h)) - exp(r*(z3+h)));

fe11(r)= (exp(r))*V11;
fe12(r)= (exp(r))*V12;
fe13(r)= (exp(r))*V13;
%--------------------------------------------------------------------------
%f11
%--------------------------------------------------------------------------
Fe11= (fe11) - Serii(@(r)fe11,0.0001,kmax); %new
F1e11= (Fe11)*(r^-0.5)*(exp(-1i*w*r))*(whittakerW(0,0,-2*1i*w*r))*exp(-r);
Saay1 = (exp(-pi*1i/4))*(subs(F1e11,r,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e11,r,-1i*t/w));
for ii=1:n
    saay_wt1 = saay_wt1 + A(ii)*subs(Saay1,t,taw(ii));
end
if11(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt1) + Serii(@(r)fe11,0,kmax)*(F1to2(kmax, w));
G1(x)= subs(if11, w, sqrt(x^2+y^2));
%-------------------------------------------------------------------------%
%f12
%--------------------------------------------------------------------------
Fe12= (fe12) - Serii(@(r)fe12,0.0001,kmax); %new
F1e12= (Fe12)*(r^-0.5)*(exp(-1i*w*r))*(whittakerW(0,0,-2*1i*w*r))*exp(-r);
Saay2 = (exp(-pi*1i/4))*(subs(F1e12,r,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e12,r,-1i*t/w));
for ii=1:n
    saay_wt2 = saay_wt2 + A(ii)*subs(Saay2,t,taw(ii));
end
if12(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt2) + Serii(@(r)fe12,0,kmax)*(F1to2(kmax, w));
G2(x)= subs(if12, w, sqrt(x^2+y^2));
%-------------------------------------------------------------------------%
%f13
%--------------------------------------------------------------------------
Fe13= (fe13) - Serii(@(r)fe13,0.0001,kmax); %new
F1e13= (Fe13)*(r^-0.5)*(exp(-1i*w*r))*(whittakerW(0,0,-2*1i*w*r))*exp(-r);
Saay3 = (exp(-pi*1i/4))*(subs(F1e13,r,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e13,r,-1i*t/w));
for ii=1:n
    saay_wt3 = saay_wt3 + A(ii)*subs(Saay3,t,taw(ii));
end
if13(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt3) + Serii(@(r)fe13,0,kmax)*(F1to2(kmax, w));
G3(x)= subs(if13, w, sqrt(x^2+y^2));
%-------------------------------------------------------------------------%
%% MOM Implimentation
x1=10* 1e-3; x2=10.5* 1e-3; x3=20.5* 1e-3;
y1=0.06* 1e-3; y2=0.17 * 1e-3; y3=0.23 * 1e-3;
N= 10; Z1=zeros(N,N); Z2=zeros(N,N);
dx1=(x1-0)/N; dy1=(y3-y2)/N; X1= 0:dx1:x1; Y1= y2:dy1:y3;
dx2=(x3-x2)/N; dy2=(y1-0)/N; X2= x1:dx2:x3; Y2= 0:dy2:y1;
%-- first eq.
for ii=1:N
   for jj=1:N
      Z1(ii,jj)= 
       
   end
end
