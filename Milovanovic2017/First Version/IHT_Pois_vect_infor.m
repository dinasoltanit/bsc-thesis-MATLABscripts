clc; clear; close all;
%% IHT of the GF
syms r w t
lex=0.01; lem=0.01; gap=0;
% fs=1e-5; % r= [fs:fs:lex lex+fs:fs:lex+gap lex+gap+fs:fs:lex+gap+lem-fs];
y=1;
epsr=2.7; epsp=1-25*1j; h=0.17; %d=0.003;
taw = [0.01689, 0.2206, 0.9232, 2.436, 5.125, 9.745];
A = [0.06166, 0.3204, 0.4197, 0.1773, 0.02055, 0.0003617];
n=6;
kmax=2*n-1; 
Nx=401; Nz=14;
%Z=zeros(Nx,Nz);
Z=sym(zeros(Nx,Nz));
%% Generating the Green's Function
%% MOM Implimentation
dx=0.05; dz=0.01; %mm
x1= [0:dx:10 10+dx:dx:20]; z1= [ -h:dz:(-h+0.06)  0:dz:0.06 ]; %mm  !!! >>>> gap is zero
%x1n=0:dx:10; z1n=0.17:dz:0.23; %mm
%-- first eq.
tic
for mx1=1:401
   for zz1=1:14 
    %for nx1=1:401
       for zp1=1:14
          if (z1(zp1) <= z1(zz1))
                %nz1 == z_prime
                %% 1
                saay_wt1=0;
                d=z1(zp1);
                % finding coeffs.
                A1= (exp(-r*d))*((epsp* cosh(r*d) * sinh(r*h) + epsr * cosh(r*h) * sinh(r*d))/((2*pi*r)*(epsr* cosh(r*h) + epsp* sinh(r*h))));
                % finding the green's function in freq domain
                V11= A1* exp(-r*(z1(zz1)-d));% for z > z_prime
                fe11= (exp(r))*V11;
                %f11
                Fe11= (fe11) - Serii(@(r)fe11,0.0001,kmax); %new
                F1e11= (Fe11)*(r^-0.5)*(exp(-1i*w*r))*(whittakerW(0,0,-2*1i*w*r))*exp(-r);
                Saay1 = (exp(-pi*1i/4))*(subs(F1e11,r,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e11,r,-1i*t/w));
                for ii=1:n
                  saay_wt1 = saay_wt1 + A(ii)*subs(Saay1,t,taw(ii));
                end
                if11=(((1i/w)/(sqrt(2*pi*w)))* saay_wt1) + Serii(@(r)fe11,0.0001,kmax)*(F1to2(kmax, w));
                % Z matrix
                Z(mx1,zz1)= subs(if11, w, sqrt(x1(mx1)^2+1));
          
          elseif ( 0 <= z1(zz1) < z1(zp1) )
                %nz1 == z_prime
                %% 2
                saay_wt2=0;
                d=z1(zp1);
                B = (exp(-r*d))*((epsp * sinh(r*h) - epsr * cosh(r*h))/((4*pi*r)*(epsr* cosh(r*h) + epsp* sinh(r*h))));
                C = (exp(-r*d))/(4*pi*r);
                V12= B * exp(-r*z1(1,zz1)) + C* exp(r*z1(1,zz1)); % for z_prime > z > 0
                fe12= (exp(r))*V12;
                %f12
                Fe12= (fe12) - Serii(@(r)fe12,0.0001,kmax); %new
                F1e12= (Fe12)*(r^-0.5)*(exp(-1i*w*r))*(whittakerW(0,0,-2*1i*w*r))*exp(-r);
                Saay2 = (exp(-pi*1i/4))*(subs(F1e12,r,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e12,r,-1i*t/w));
                for ii=1:n
                  saay_wt2 = saay_wt2 + A(ii)*( subs(Saay2,t,taw(ii)) );
                end
                if12=(((1i/w)/(sqrt(2*pi*w)))* saay_wt2) + Serii(@(r)fe12,0.0001,kmax)*(F1to2(kmax, w));
                Z(mx1,zz1)= subs(if12, w, (sqrt(x1(mx1)^2+1)) );
          
          elseif ( -h < z1(zz1) < 0 )
                %nz1 == z_prime
                %% 3
                saay_wt3=0;
                d=z1(zp1);
                D = (-epsp * exp(-r*d))/((4*pi*r)*(epsr* cosh(r*h) + epsp* sinh(r*h)));
                V13= D * (exp(-r*(z1(zz1)+h)) - exp(r*(z1(zz1)+h))); % for 0 > z > -h
                fe13= (exp(r))*V13;
                %f13
                Fe13= (fe13) - Serii(@(r)fe13,0.0001,kmax); %new
                F1e13= (Fe13)*(r^-0.5)*(exp(-1i*w*r))*(whittakerW(0,0,-2*1i*w*r))*exp(-r);
                Saay3= (exp(-pi*1i/4))*(subs(F1e13,r,1i*t/w)) - (exp(pi*1i/4))*(subs(F1e13,r,-1i*t/w));
                for ii=1:n
                  saay_wt3 = saay_wt3 + A(ii)*subs(Saay3,t,taw(ii));
                end
                if13=(((1i/w)/(sqrt(2*pi*w)))* saay_wt3) + Serii(@(r)fe13,0.0001,kmax)*(F1to2(kmax, w));
                Z(mx1,zz1)= subs(if13, w, sqrt(x1(mx1)^2+1));
          end
       end
    %end
   end
end
toc






