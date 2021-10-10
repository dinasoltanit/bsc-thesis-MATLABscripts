clc;clear;
%-----------------------Parameters--------------------------------------------------------------------------------------------------
syms x z w t r1 r2
eps1=1.2;
eps2=1.5;
h1=-0.2;
h2=-0.4;
z_prime=5;
%-----------------------Defining the solution: v_tilda------------------------------------------------------------------------------
fA =(eps1*exp(2*x*(h1 + h2)) + eps2*exp(2*x*(h1 + h2)) - eps1*exp(2*x*(h1 + z_prime)) + eps1*exp(2*x*(h2 + z_prime)) - eps2*exp(2*x*(h1 + z_prime)) - eps2*exp(2*x*(h2 + z_prime)) - eps1*exp(2*h1*x) + eps1*exp(2*h2*x) - eps2*exp(2*h1*x) - eps1*exp(4*h1*x) - eps2*exp(2*h2*x) + eps2*exp(4*h1*x) - eps1^2*exp(2*x*(h1 + h2)) - eps1*exp(2*x*(2*h1 + z_prime)) + eps2*exp(2*x*(2*h1 + z_prime)) + eps1^2*exp(2*x*(h1 + z_prime)) - eps1^2*exp(2*x*(h2 + z_prime)) + eps1*exp(2*x*(h1 + h2 + z_prime)) + eps2*exp(2*x*(h1 + h2 + z_prime)) - eps1^2*exp(2*h1*x) + eps1^2*exp(2*h2*x) + eps1^2*exp(4*h1*x) - eps1^2*exp(2*x*(2*h1 + z_prime)) + eps1^2*exp(2*x*(h1 + h2 + z_prime)) + eps1*eps2*exp(2*x*(h1 + h2 + z_prime)) - eps1*eps2*exp(2*x*(h1 + h2)) + eps1*eps2*exp(2*x*(h1 + z_prime)) + eps1*eps2*exp(2*x*(h2 + z_prime)) - eps1*eps2*exp(2*h1*x) - eps1*eps2*exp(2*h2*x) - eps1*eps2*exp(4*h1*x) + eps1*eps2*exp(2*x*(2*h1 + z_prime)))/(4*x*pi*(eps1*exp(2*x*(h2 + z_prime)) - eps1*exp(2*x*(h1 + z_prime)) - eps2*exp(2*x*(h1 + z_prime)) - eps2*exp(2*x*(h2 + z_prime)) - eps1*exp(2*x*(2*h1 + z_prime)) + eps2*exp(2*x*(2*h1 + z_prime)) + eps1^2*exp(2*x*(h1 + z_prime)) - eps1^2*exp(2*x*(h2 + z_prime)) + eps1*exp(2*x*(h1 + h2 + z_prime)) + eps2*exp(2*x*(h1 + h2 + z_prime)) - eps1^2*exp(2*x*(2*h1 + z_prime)) + eps1^2*exp(2*x*(h1 + h2 + z_prime)) + eps1*eps2*exp(2*x*(h1 + h2 + z_prime)) + eps1*eps2*exp(2*x*(h1 + z_prime)) + eps1*eps2*exp(2*x*(h2 + z_prime)) + eps1*eps2*exp(2*x*(2*h1 + z_prime))));
fB1 = -(eps1*exp(2*h1*x) - eps2*exp(2*x*(h1 + h2)) - eps1*exp(2*x*(h1 + h2)) - eps1*exp(2*h2*x) + eps2*exp(2*h1*x) + eps1*exp(4*h1*x) + eps2*exp(2*h2*x) - eps2*exp(4*h1*x) + eps1^2*exp(2*x*(h1 + h2)) + eps1^2*exp(2*h1*x) - eps1^2*exp(2*h2*x) - eps1^2*exp(4*h1*x) + eps1*eps2*exp(2*x*(h1 + h2)) + eps1*eps2*exp(2*h1*x) + eps1*eps2*exp(2*h2*x) + eps1*eps2*exp(4*h1*x))/(4*x*pi*(eps1*exp(x*(2*h1 + 2*h2 + z_prime)) + eps2*exp(x*(2*h1 + 2*h2 + z_prime)) - eps1*exp(x*(2*h1 + z_prime)) + eps1*exp(x*(2*h2 + z_prime)) - eps2*exp(x*(2*h1 + z_prime)) - eps1*exp(x*(4*h1 + z_prime)) - eps2*exp(x*(2*h2 + z_prime)) + eps2*exp(x*(4*h1 + z_prime)) + eps1^2*exp(x*(2*h1 + 2*h2 + z_prime)) + eps1^2*exp(x*(2*h1 + z_prime)) - eps1^2*exp(x*(2*h2 + z_prime)) - eps1^2*exp(x*(4*h1 + z_prime)) + eps1*eps2*exp(x*(2*h1 + 2*h2 + z_prime)) + eps1*eps2*exp(x*(2*h1 + z_prime)) + eps1*eps2*exp(x*(2*h2 + z_prime)) + eps1*eps2*exp(x*(4*h1 + z_prime))));
fB2 = exp(-x*z_prime)/(4*x*pi); 
fC1 =-(eps1*exp(3*h1*x) + eps2*exp(3*h1*x) - eps1*exp(x*(h1 + 2*h2)) + eps2*exp(x*(h1 + 2*h2)))/(2*x*pi*(eps1*exp(x*(2*h1 + 2*h2 + z_prime)) + eps2*exp(x*(2*h1 + 2*h2 + z_prime)) - eps1*exp(x*(2*h1 + z_prime)) + eps1*exp(x*(2*h2 + z_prime)) - eps2*exp(x*(2*h1 + z_prime)) - eps1*exp(x*(4*h1 + z_prime)) - eps2*exp(x*(2*h2 + z_prime)) + eps2*exp(x*(4*h1 + z_prime)) + eps1^2*exp(x*(2*h1 + 2*h2 + z_prime)) + eps1^2*exp(x*(2*h1 + z_prime)) - eps1^2*exp(x*(2*h2 + z_prime)) - eps1^2*exp(x*(4*h1 + z_prime)) + eps1*eps2*exp(x*(2*h1 + 2*h2 + z_prime)) + eps1*eps2*exp(x*(2*h1 + z_prime)) + eps1*eps2*exp(x*(2*h2 + z_prime)) + eps1*eps2*exp(x*(4*h1 + z_prime))));
fC2 =(eps2*exp(3*h1*x) - eps1*exp(3*h1*x) + eps1*exp(x*(h1 + 2*h2)) + eps2*exp(x*(h1 + 2*h2)))/(2*x*pi*(eps1*exp(x*(2*h1 + 2*h2 + z_prime)) + eps2*exp(x*(2*h1 + 2*h2 + z_prime)) - eps1*exp(x*(2*h1 + z_prime)) + eps1*exp(x*(2*h2 + z_prime)) - eps2*exp(x*(2*h1 + z_prime)) - eps1*exp(x*(4*h1 + z_prime)) - eps2*exp(x*(2*h2 + z_prime)) + eps2*exp(x*(4*h1 + z_prime)) + eps1^2*exp(x*(2*h1 + 2*h2 + z_prime)) + eps1^2*exp(x*(2*h1 + z_prime)) - eps1^2*exp(x*(2*h2 + z_prime)) - eps1^2*exp(x*(4*h1 + z_prime)) + eps1*eps2*exp(x*(2*h1 + 2*h2 + z_prime)) + eps1*eps2*exp(x*(2*h1 + z_prime)) + eps1*eps2*exp(x*(2*h2 + z_prime)) + eps1*eps2*exp(x*(4*h1 + z_prime)))); 
fD = -(eps1*exp(x*(2*h1 + h2)))/(x*pi*(eps1*exp(x*(2*h1 + 2*h2 + z_prime)) + eps2*exp(x*(2*h1 + 2*h2 + z_prime)) - eps1*exp(x*(2*h1 + z_prime)) + eps1*exp(x*(2*h2 + z_prime)) - eps2*exp(x*(2*h1 + z_prime)) - eps1*exp(x*(4*h1 + z_prime)) - eps2*exp(x*(2*h2 + z_prime)) + eps2*exp(x*(4*h1 + z_prime)) + eps1^2*exp(x*(2*h1 + 2*h2 + z_prime)) + eps1^2*exp(x*(2*h1 + z_prime)) - eps1^2*exp(x*(2*h2 + z_prime)) - eps1^2*exp(x*(4*h1 + z_prime)) + eps1*eps2*exp(x*(2*h1 + 2*h2 + z_prime)) + eps1*eps2*exp(x*(2*h1 + z_prime)) + eps1*eps2*exp(x*(2*h2 + z_prime)) + eps1*eps2*exp(x*(4*h1 + z_prime))));
%--------------------------------------------------------------------------
taw = [0.01689, 0.2206, 0.9232, 2.436, 5.125, 9.745];
A = [0.06166, 0.3204, 0.4197, 0.1773, 0.02055, 0.0003617];
n=6;
kmax=2*n-1;
saay_wt1=0;
saay_wt2=0;
saay_wt3=0;
saay_wt4=0;

w2=(1.1:0.1:6);

z1=(5:0.2:9.8);
z2=(0:0.2:4.8);
z3=(-0.2:0.008:-0.008);
z4=(-0.4:0.008:-0.208);

if11Vect1=zeros(1,25);
if12Vect1=zeros(1,25);
if13Vect1=zeros(1,25);
if14Vect1=zeros(1,25);
%--------------------------------------------------------------------------
V11(x)= x*(fA*exp(-k*z));
V12(x)= x*(fB1*exp(-k*(z-z_prime)) + fB2*exp(k*(z-z_prime)));
V13(x)= x*(fC1*exp(-k*(z+h1)) + fC2*exp(k*(z+h1)));
V14(x)= x*(fD*(exp(-k*(z+h2)) - exp(k*(z+h2))));
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fe11z(x)= (exp(x))*V11;
fe12z(x)= (exp(x))*V12;
fe13z(x)= (exp(x))*V13;
fe14z(x)= (exp(x))*V14;
%--------------------------------------------------------------------------
%f11
%--------------------------------------------------------------------------
for jj=1:25
    fe11=subs(fe11z,z,z1(jj));
    for kk=1:25
        Fe11= (fe11) - Serii(@(x)fe11,0.0001,kmax); %new
        F1e11= (Fe11)*(x^-0.5)*(exp(-1i*w*x))*(whittakerW(0,0,-2*1i*w*x))*exp(-x);
        Saay1 = (exp(-pi*1i/4))*(subs(F1e11,x,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e11,x,-1i*t/w));
            for ii=1:n
                saay_wt1 = saay_wt1 + A(ii)*subs(Saay1,t,taw(ii));
            end
        if11(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt1) + Serii(@(x)fe11,0,kmax)*(F1to2(kmax, w));
        %I1=simplify(subs(if11,w,sqrt(r1^2+r2^2)));
        if11Vect1(jj,kk)=(subs(if11,w,w2(kk)));
    end
end
%-------------------------------------------------------------------------
%f12
%--------------------------------------------------------------------------
for jj=1:25
    fe12=subs(fe12z,z,z2(jj));
    for kk=1:25
    Fe12= (fe12) - Serii(@(x)fe12,0.0001,kmax); %new
    F1e12= (Fe12)*(x^-0.5)*(exp(-1i*w*x))*(whittakerW(0,0,-2*1i*w*x))*exp(-x);
    Saay2 = (exp(-pi*1i/4))*(subs(F1e12,x,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e12,x,-1i*t/w));
        for ii=1:n
            saay_wt2 = saay_wt2 + A(ii)*subs(Saay2,t,taw(ii));
        end
    if12(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt2) + Serii(@(x)fe12,0,kmax)*(F1to2(kmax, w));
    %I2=simplify(subs(if11,w,sqrt(r1^2+r2^2)));
    if12Vect1(jj,kk)=(subs(if12,w,w2(kk)));
    end
end
%-------------------------------------------------------------------------%
%f13
%--------------------------------------------------------------------------
for jj=1:25
    fe13=subs(fe13z,z,z3(jj));
    for kk=1:25
    Fe13= (fe13) - Serii(@(x)fe13,0.0001,kmax); %new
    F1e13= (Fe13)*(x^-0.5)*(exp(-1i*w*x))*(whittakerW(0,0,-2*1i*w*x))*exp(-x);
    Saay3 = (exp(-pi*1i/4))*(subs(F1e13,x,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e13,x,-1i*t/w));
        for ii=1:n
            saay_wt3 = saay_wt3 + A(ii)*subs(Saay3,t,taw(ii));
        end
    if13(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt3) + Serii(@(x)fe13,0,kmax)*(F1to2(kmax, w));
    %I3=simplify(subs(if13,w,sqrt(r1^2+r2^2)));
    if13Vect1(jj,kk)=(subs(if13,w,w2(kk)));
    end
end
%-------------------------------------------------------------------------%
%f14
%--------------------------------------------------------------------------
for jj=1:25
    fe14=subs(fe14z,z,z4(jj));
    for kk=1:25
    Fe14= (fe14) - Serii(@(x)fe14,0.0001,kmax); %new
    F1e14= (Fe14)*(x^-0.5)*(exp(-1i*w*x))*(whittakerW(0,0,-2*1i*w*x))*exp(-x);
    Saay4 = (exp(-pi*1i/4))*(subs(F1e14,x,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e14,x,-1i*t/w));
        for ii=1:n
            saay_wt4 = saay_wt4 + A(ii)*subs(Saay4,t,taw(ii));
        end
    if14(w)=(((1i/w)/(sqrt(2*pi*w)))* saay_wt4) + Serii(@(x)fe14,0,kmax)*(F1to2(kmax, w));
    %I4=simplify(subs(if14,w,sqrt(r1^2+r2^2)));
    if14Vect1(jj,kk)=(subs(if14,w,w2(kk)));
    end
end
GF = piecewise(z>z_prime,I1, (0<z)&(z<z_prime), I2, -h1<z<0, I3, -h2<z<-h1, I4);
%---------------------------------Plots------------------------------------
subplot(1,4,1)
surf(w2,z1,abs(if11Vect1))
grid on;

subplot(1,4,2)
surf(w2,z2,abs(if12Vect1))
grid on;

subplot(1,4,3)
surf(w2,z3,abs(if13Vect1))
grid on;

subplot(1,4,4)
surf(w2,z4,abs(if14Vect1))