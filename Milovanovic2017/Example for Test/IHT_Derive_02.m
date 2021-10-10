clc;
clear;
syms x w t
z1=7.5;
z2=2.5;
z3=0.1;
z4=0.3;
taw = [0.01689, 0.2206, 0.9232, 2.436, 5.125, 9.745];
A = [0.06166, 0.3204, 0.4197, 0.1773, 0.02055, 0.0003617];
n=6;
kmax=2*n-1;
saay_wt1=0;
saay_wt2=0;
saay_wt3=0;
saay_wt4=0;
%--------------------------------------------------------------------------
V11(x)= -x*((exp(-x*(z1 - 5))*((297*exp((2*x)/5))/50 + (18*exp((4*x)/5))/25 + (27*exp((6*x)/5))/50 - (27*exp((52*x)/5))/50 - (18*exp((54*x)/5))/25 - (297*exp((56*x)/5))/50))/(4*x*pi*((27*exp((52*x)/5))/50 + (18*exp((54*x)/5))/25 + (297*exp((56*x)/5))/50)));
V12(x)= x*((exp(-5*x)*exp(x*z2))/(4*x*pi) - (exp(-x*z2)*((297*exp((2*x)/5))/50 + (18*exp((4*x)/5))/25 + (27*exp((6*x)/5))/50))/(4*x*pi*((27*exp((27*x)/5))/50 + (18*exp((29*x)/5))/25 + (297*exp((31*x)/5))/50)));
V13(x) = x*((exp(x*(z3 + 1/5))*((3*exp((3*x)/5))/10 + (27*exp(x))/10))/(2*x*pi*((27*exp((27*x)/5))/50 + (18*exp((29*x)/5))/25 + (297*exp((31*x)/5))/50)) - (exp(-x*(z3 + 1/5))*((27*exp((3*x)/5))/10 + (3*exp(x))/10))/(2*x*pi*((27*exp((27*x)/5))/50 + (18*exp((29*x)/5))/25 + (297*exp((31*x)/5))/50)));
V14(x) = -x*((6*exp((4*x)/5)*(exp(x*(z4 + 2/5)) + exp(-x*(z4 + 2/5))))/(5*x*pi*((27*exp((27*x)/5))/50 + (18*exp((29*x)/5))/25 + (297*exp((31*x)/5))/50)));
%--------------------------------------------------------------------------
fe11(x)= (exp(x))*V11;
fe12(x)= (exp(x))*V12;
fe13(x)= (exp(x))*V13;
fe14(x)= (exp(x))*V14;
%--------------------------------------------------------------------------
%f11
%--------------------------------------------------------------------------
Fe11= (fe11) - Serii(@(x)fe11,0.0001,kmax); %new
F1e11= (Fe11)*(x^-0.5)*(exp(-1i*w*x))*(whittakerW(0,0,-2*1i*w*x))*exp(-x);
Saay1 = (exp(-pi*1i/4))*(subs(F1e11,x,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e11,x,-1i*t/w));
for ii=1:n
    saay_wt1 = saay_wt1 + A(ii)*subs(Saay1,t,taw(ii));
end
if11(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt1) + Serii(@(x)fe11,0,kmax)*(F1to2(kmax, w));
%I1=subs(if11,w,sqrt(r1,r2));
%-------------------------------------------------------------------------%
%f12
%--------------------------------------------------------------------------
Fe12= (fe12) - Serii(@(x)fe12,0.0001,kmax); %new
F1e12= (Fe12)*(x^-0.5)*(exp(-1i*w*x))*(whittakerW(0,0,-2*1i*w*x))*exp(-x);
Saay2 = (exp(-pi*1i/4))*(subs(F1e12,x,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e12,x,-1i*t/w));
for ii=1:n
    saay_wt2 = saay_wt2 + A(ii)*subs(Saay2,t,taw(ii));
end
if12(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt2) + Serii(@(x)fe12,0,kmax)*(F1to2(kmax, w));
%I2=subs(if12,w,sqrt(r1,r2));
%-------------------------------------------------------------------------%
%f13
%--------------------------------------------------------------------------
Fe13= (fe13) - Serii(@(x)fe13,0.0001,kmax); %new
F1e13= (Fe13)*(x^-0.5)*(exp(-1i*w*x))*(whittakerW(0,0,-2*1i*w*x))*exp(-x);
Saay3 = (exp(-pi*1i/4))*(subs(F1e13,x,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e13,x,-1i*t/w));
for ii=1:n
    saay_wt3 = saay_wt3 + A(ii)*subs(Saay3,t,taw(ii));
end
if13(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt3) + Serii(@(x)fe13,0,kmax)*(F1to2(kmax, w));
%I3=subs(if13,w,sqrt(r1,r2));
%-------------------------------------------------------------------------%
%f14
%--------------------------------------------------------------------------
Fe14= (fe14) - Serii(@(x)fe14,0.0001,kmax); %new
F1e14= (Fe14)*(x^-0.5)*(exp(-1i*w*x))*(whittakerW(0,0,-2*1i*w*x))*exp(-x);
Saay4 = (exp(-pi*1i/4))*(subs(F1e14,x,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e14,x,-1i*t/w));
for ii=1:n
    saay_wt4 = saay_wt4 + A(ii)*subs(Saay4,t,taw(ii));
end
if14(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt4) + Serii(@(x)fe14,0,kmax)*(F1to2(kmax, w));
%I4=subs(if14,w,sqrt(r1,r2));
%-------------------------------------------------------------------------%
%w1=(-5:-1);
w2=(1.1:0.1:6);
if11Vect1=zeros(1,50);
if11Vect2=zeros(1,50);
if12Vect1=zeros(1,50);
if12Vect2=zeros(1,50);
if13Vect1=zeros(1,50);
if13Vect2=zeros(1,50);
if14Vect1=zeros(1,50);
if14Vect2=zeros(1,50);
% it is finding the value of the v function at specific freqencies
for ii=1:50
    if11Vect1(1,ii)=(subs(if11,w,w2(ii)));   
end

% for ii=1:50
%     if11Vect2(1,ii)=(subs(if11,w,w1(ii)));   
% end

for ii=1:50
    if12Vect1(1,ii)=(subs(if12,w,w2(ii)));   
end

% for ii=1:50
%     if12Vect2(1,ii)=(subs(if12,w,w1(ii)));   
% end

for ii=1:50
    if13Vect1(1,ii)=(subs(if13,w,w2(ii)));   
end

% for ii=1:50
%     if13Vect2(1,ii)=(subs(if13,w,w1(ii)));   
% end

for ii=1:50
    if14Vect1(1,ii)=(subs(if14,w,w2(ii)));   
end

% for ii=1:50
%     if14Vect2(1,ii)=(subs(if14,w,w1(ii)));   
% end

% subplot(4,2,1)
% plot(w1,abs(if11Vect2),'g')
% grid on;
% title('z=7.5 , rho=-5:-1');
% xlabel('rho'); ylabel('Potential');
%--------------------------------------------------------------------------
subplot(4,1,1)
plot(w2,abs(if11Vect1),'g')
grid on;
title('z=7.5 , rho=1.1:0.1:6');
xlabel('rho'); ylabel('Potential');

% subplot(4,2,3)
% plot(w1,abs(if12Vect2),'g')
% grid on;
% title('z=2.5 , rho=-5:-1');
% xlabel('rho'); ylabel('Potential');

subplot(4,1,2)
plot(w2,abs(if12Vect1),'g')
grid on;
title('z=2.5 , rho=1.1:0.1:6');
xlabel('rho'); ylabel('Potential');

% subplot(4,2,5)
% plot(w1,abs(if13Vect2),'g')
% grid on;
% title('z=-0.1 , rho=-5:-1');
% xlabel('Frequency (Hz)'); ylabel('Amplitude');

subplot(4,1,3)
plot(w2,abs(if13Vect1),'g')
grid on;
title('z=-0.1 , rho=1.1:0.1:6');
xlabel('rho'); ylabel('Potential');

% subplot(4,2,7)
% plot(w1,abs(if14Vect2),'g')
% grid on;
% title('z=-0.3 , rho=-5:-1');
% xlabel('rho'); ylabel('Potential');

subplot(4,1,4)
plot(w2,abs(if14Vect1),'g')
grid on;
title('z=-0.3 , rho=1.1:0.1:6');
xlabel('rho'); ylabel('Potential');