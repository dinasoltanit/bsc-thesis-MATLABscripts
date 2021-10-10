clc;
clear;

%-----------------------------General--------------------------------------
ref=0.05; %mm

%-----------------------------Dimensions-----------------------------------
h1=0.1; %mm
l1=0.5;
%--
h2=0.1; %mm
l2=3;
%--
a=1.5;
b=3; %mm %h3,l3
%--
h4=1; %mm
l4=0.25;
%--
h5=0.1; %mm
l5=0.75;
%--
l6=21.5;
%--
h8=1.1; %mm
l8=13;
%--
h9=1.1; %mm
l9=4.75;

%--------------------- dielectric-to-dielectric interfaces ----------------
da1=l8/ref;
da2=l9/ref;
da3=a/ref;
da4=(sqrt(a^2+b^2))/ref;
da5=l4/ref;
%--
dd1=b/ref;
dd2=h4/ref;
dd3=h4/ref;
dd4=l5/ref;
%--------------------- dielectric-to-conductor interfaces -----------------
dc1=l1/ref;
dc2=l2/ref;
dc3=h2/ref;
dc4=h2/ref;
dc5=h1/ref;
%--
ca1=l1/ref;
ca2=h1/ref;
%-------------------------------- ground ----------------------------------
gnd=l6/ref;
%-------------------------------- ground ----------------------------------
Ns= dc1+dc2+dc3+dc4+dc5+ca1+ca2; %without gnd-to-conductor interface
Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;