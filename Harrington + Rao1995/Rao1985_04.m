clc;
clear;
% if * then not-checked if ** then checked.
%syms Vp Vn ep
%-----------------------------General-------------------------------------**
ref=0.05;   nh=11.1/ref;    nl=21.5/ref;    dx=ref; dy=ref; Vp=5000;    Vn=0;
eps0=8.85*10^(-12); eps1=1.5;   eps2=4; eps3=1+1j;  %eps of plasma which is a function
deltat=ref; K=1;
%-----------------------------Dimensions----------------------------------*
h1=0.1; %mm
l1=0.5;
%--
h2=0.1; %mm
l2=3;
%--
a=1.5;
h3=a-h1;
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
da1=l8/ref; da2=l9/ref; da3=a/ref;  da4=a/ref;  da5=l4/ref;
%--
dd1=b/ref;  dd2=h4/ref; dd3=h4/ref; dd4=l5/ref;
%--------------------- dielectric-to-conductor interfaces -----------------
cd1=l1/ref; cd2=l2/ref; cd3=h2/ref; cd4=h2/ref; cd5=h1/ref;
%--
ca1=l1/ref; ca2=h1/ref;
%-------------------------------- ground ----------------------------------
gnd=l6/ref;
%-------------------------------- ground ----------------------------------
Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2; %without gnd-to-conductor interface
Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
N= Ns + Nd;
%-------------------------------- Coordinate Matrix ----------------------**
%[X,Y]=meshgrid(0:0.05:21.5,11.1:-0.05:0);
% Z=sym(zeros(N,N));
% V=sym(zeros(N,1));
Z=zeros(N,N);
V=zeros(N,1);
%sigma=sym(zeros(N,1));
sigma=zeros(N,1);
%-------------------------------- Solution Procedure-----------------------
% rs=zeros(Ns,3);
% rn=zeros(Nd,3);
r=zeros(N,3);
rhoi12=zeros(N,2); %r3
lhat=zeros(N,2);
nhat=zeros(N,2);
xi12=zeros(N,1);
yi12=zeros(N,1);
%-------------------------------- Deriving The Z Matrix -------------------
%-- j=1,2,...,Ns
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%% ----------------------------- Deriving rs Matrix ----------------------**
% cd1
for j=1:cd1
    %rs=[j,xj,yj]
    r(j,1)=j;
    r(j,2)=13+(j-1)*0.05; r(j,3)=1.1;
    %% i1i2
    if j==cd1
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==cd1
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
% cd2
for j=1+cd1:cd2+cd1
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13.75+(j-1)*0.05; r(j,3)=0.1;
    %% i1i2
    if j==cd2+cd1
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==cd2+cd1
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%cd3
for j=1+cd1+cd2:cd3+cd1+cd2
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=16.75; r(j,3)=0+(j-1)*0.05;
    %% i1i2
    if j==cd3+cd1+cd2
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==cd3+cd1+cd2
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%cd4
for j=1+cd1+cd2+cd3:cd4+cd1+cd2+cd3
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13.75; r(j,3)=0+(j-1)*0.05;
    %% i1i2
    if j==cd4+cd1+cd2+cd3
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==cd4+cd1+cd2+cd3
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%cd5
for j=1+cd1+cd2+cd3+cd4:cd5+cd1+cd2+cd3+cd4
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13.5; r(j,3)=1.1+(j-1)*0.05;
    %% i1i2
    if j==cd5+cd1+cd2+cd3+cd4
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==cd5+cd1+cd2+cd3+cd4
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%ca1
for j=1+cd1+cd2+cd3+cd4+cd5:ca1+cd1+cd2+cd3+cd4+cd5
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13+(j-1)*0.05; r(j,3)=1.2;
    %% i1i2
    if j==ca1+cd1+cd2+cd3+cd4+cd5
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==ca1+cd1+cd2+cd3+cd4+cd5
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%ca2
for j=1+cd1+cd2+cd3+cd4+cd5+ca1:ca2+cd1+cd2+cd3+cd4+cd5+ca1
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13; r(j,3)=1.1+(j-1)*0.05;
    %% i1i2
    if j==ca2+cd1+cd2+cd3+cd4+cd5+ca1
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==ca2+cd1+cd2+cd3+cd4+cd5+ca1
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%% ----------------------------- Deriving rd Matrix ----------------------**
%r=[i,xi,yi]
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da1
for j=1+Ns:da1+Ns
   r(j,1)=j; r(j,2)=0+(j-1)*0.05; r(j,3)=1.1;
   %% i1i2
    if j==da1+Ns
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==da1+Ns
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da2
for j=1+Ns+da1:da2+Ns+da1
   r(j,1)=j; r(j,2)=16.75+(j-1)*0.05; r(j,3)=1.1;
   %% i1i2
    if j==da2+Ns+da1
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==da2+Ns+da1
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da3
for j=1+Ns+da1+da2:da3+Ns+da1+da2
   r(j,1)=j; r(j,2)=13.5; r(j,3)=1.2+(j-1)*0.05;
   %% i1i2
    if j==da3+Ns+da1+da2
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==da3+Ns+da1+da2
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da4
for i=1+Ns+da1+da2+da3:da4+Ns+da1+da2+da3
    r(j,1)=j; r(j,2)=13.5+(j-1)*(b/da4); r(j,3)=1.2+(j-1)*(a/da4);
    %% i1i2
    if j==da4+Ns+da1+da2+da3
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==da4+Ns+da1+da2+da3
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da5
for j=1+Ns+da1+da2+da3+da4:da5+Ns+da1+da2+da3+da4
   r(j,1)=j; r(j,2)=16.5+(j-1)*0.05; r(j,3)=1.1;
   %% i1i2
    if j==da5+Ns+da1+da2+da3+da4
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==da5+Ns+da1+da2+da3+da4
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd1
for j=1+Ns+da1+da2+da3+da4+da5:dd1+Ns+da1+da2+da3+da4+da5
   r(j,1)=j; r(j,2)=13.5+(j-1)*0.05; r(j,3)=1.1;
   %% i1i2
    if j==dd1+Ns+da1+da2+da3+da4+da5
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==dd1+Ns+da1+da2+da3+da4+da5
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd2
for j=1+Ns+da1+da2+da3+da4+da5+dd1:dd2+Ns+da1+da2+da3+da4+da5+dd1
   r(j,1)=j; r(j,2)=16.75; r(j,3)=0.1+(j-1)*0.05;
   %% i1i2
    if j==dd2+Ns+da1+da2+da3+da4+da5+dd1
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==dd2+Ns+da1+da2+da3+da4+da5+dd1
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd3
for j=1+Ns+da1+da2+da3+da4+da5+dd1+dd2:dd3+Ns+da1+da2+da3+da4+da5+dd1+dd2
   r(j,1)=j; r(j,2)=13; r(j,3)=0.1+(j-1)*0.05;
   %% i1i2
    if j==dd3+Ns+da1+da2+da3+da4+da5+dd1+dd2
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==dd3+Ns+da1+da2+da3+da4+da5+dd1+dd2
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd4
for j=1+Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3:dd4+Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
   r(j,1)=j; r(j,2)=13+(j-1)*0.05; r(j,3)=0.1;
   %% i1i2
    if j==dd4+Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
        xi12(j,1)=(r(j,2) + 0.05)-r(j,2);
    else
        xi12(j,1)=(r(j+1,2))-r(j,2);
    end
    
    if j==dd4+Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
        yi12(j,1)=(r(j,3) + 0.05)-r(j,3);
    else
        yi12(j,1)=(r(j+1,3) + 0.05)-r(j,3);
    end
        rhoi12(j,1)=xi12(j,1); %x3
        rhoi12(j,2)=yi12(j,1); %y3
    %%  l_hat
    lhat(j,1)= xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    lhat(j,2)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    %% n_hat
    nhat(j,1)= yi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
    nhat(j,2)=-xi12(j,1)/(sqrt((rhoi12(j,1))^2+(rhoi12(j,2))^2));
end
%% ---------------------------------Vectors--------------------------------
xji1=zeros(N,N);
yji1=zeros(N,N);
xji2=zeros(N,N);
yji2=zeros(N,N);
% xi1i2=zeros(N,N);
% yi1i2=zeros(N,N);
rhoji1=zeros(N,N,2);  %r1
rhoji2=zeros(N,N,2);  %r2
% rhoi1i2=zeros(N,N,2); %r3
% lhat=zeros(N,N,2);
% nhat=zeros(N,N,2);
Rj1=zeros(N,N);   %--1
Rj1c=zeros(N,N);
Rj1l=zeros(N,N);
arg_j1=zeros(N,N);
Rj2=zeros(N,N);   %--2
Rj2c=zeros(N,N);
Rj2l=zeros(N,N);
arg_j2=zeros(N,N);
R12=zeros(N,N);   %--3
%R7=zeros(N,N);   %--4
S=zeros(N,N);     %--5
delta=zeros(N,N); %--6
R1=zeros(N,N);    %--7
RN=zeros(N,N);    %--8
RL=zeros(N,N);    %--9
L1=zeros(N,N,2);  %--10
R18= zeros(N,N);  %--11
I=zeros(Nd,N);
for j=1:N
    for i=1:N
        %ji1
        xji1(j,i)=r(j,2)-r(i,2);
        yji1(j,i)=r(j,3)-r(i,3);
        %ji2
%         if i==N
%             xji2(j,i)=r(j,2)-r(i,2);
%             yji2(j,i)=r(j,3)-r(i,3);
%         else
%             xji2(j,i)=r(j,2)-(r(i,2) + 0.05);
%             yji2(j,i)=r(j,3)-(r(i,3) + 0.05);
%         end
        %ji2
            xji2(j,i)=r(j,2)-(r(i,2) + 0.05);
            yji2(j,i)=r(j,3)-(r(i,3) + 0.05);
%         %i1i2
%         if i==N
%             xi1i2(j,i)=r(i,2)-r(i,2);
%             yi1i2(j,i)=r(i,3)-r(i,3);    
%         else
%             xi1i2(j,i)=(r(i,2) + 0.05)-r(i,2);
%             yi1i2(j,i)=(r(i,3) + 0.05)-r(i,3);
%         end
%         %i1i2
%             xi1i2(j,i)=(r(i,2) + 0.05)-r(i,2);
%             yi1i2(j,i)=(r(i,3) + 0.05)-r(i,3);
        % r1
        rhoji1(j,i,1)=xji1(j,i); %x1
        rhoji1(j,i,2)=yji1(j,i); %y1
        %r2
        rhoji2(j,i,1)=xji2(j,i); %x2
        rhoji2(j,i,2)=yji2(j,i); %y2
        %r3
%         rhoi12(j,i,1)=xi12(j,i); %x3
%         rhoi12(j,i,2)=yi12(j,i); %y3
%         % l_hat
%         lhat(j,i,1)=xi1i2(j,i)/(sqrt((rhoi1i2(j,i,1))^2+(rhoi1i2(j,i,2))^2));
%         lhat(j,i,2)=yi1i2(j,i)/(sqrt((rhoi1i2(j,i,1))^2+(rhoi1i2(j,i,2))^2));
%         % n_hat
%         nhat(j,i,1)= yi1i2(j,i)/(sqrt((rhoi1i2(j,i,1))^2+(rhoi1i2(j,i,2))^2));
%         nhat(j,i,2)=-xi1i2(j,i)/(sqrt((rhoi1i2(j,i,1))^2+(rhoi1i2(j,i,2))^2));
        %% --1 Rj1
        Rj1c(j,i)=((lhat(i,1)*xji1(j,i))+(lhat(i,2)*yji1(j,i)));
        arg_j1(j,i)=(sqrt((rhoji1(j,i,1))^2+(rhoji1(j,i,2))^2));
        Rj1l(j,i)=log(arg_j1(j,i));
        if (arg_j1(j,i)<1.0e-4)
            if (Rj1c(j,i)<1.0e-4)
                Rj1(j,i)=0;
            end
        else
           %Rj1(j,i)=log((sqrt((rhoji1(j,i,1))^2+(rhoji1(j,i,2))^2)))*((lhat(j,i,1)*xji1(j,i))+(lhat(j,i,2)*yji1(j,i)));
            Rj1(j,i)= Rj1l(j,i) * Rj1c(j,i);
        end
        %% --2 Rj2
        Rj2c(j,i)=((lhat(i,1)*xji2(j,i))+(lhat(i,2)*yji2(j,i)));
        arg_j2(j,i)=(sqrt((rhoji2(j,i,1))^2+(rhoji2(j,i,2))^2));
        Rj1l(j,i)=log(arg_j2(j,i));
        if (arg_j2(j,i)<1.0e-4)
            if (Rj2c(j,i)<1.0e-4)
                Rj2(j,i)=0;
            end
        else
            Rj2(j,i)= Rj2l(j,i) * Rj2c(j,i);
           %Rj2(j,i)=log((sqrt((rhoji2(j,i,1))^2+(rhoji2(j,i,2))^2)))*((lhat(j,i,1)*xji2(j,i))+(lhat(j,i,2)*yji2(j,i)));
        end
        %% --3 R12
        R12(j,i)= (sqrt((rhoi12(i,1))^2+(rhoi12(i,2))^2));
        %--4 R7
        %R7(j,i)= sqrt(((rhoji1(i,1))^2+(rhoji1(i,2))^2)-(((lhat(i,1)*xji1(j,i))+(lhat(i,2)*yji1(j,i)))^2));
        %--5 S
        S(j,i)= ((xji1(j,i))*(xji2(j,i)))+((yji1(j,i))*(yji2(j,i)));
        %--6 delta
        delta(j,i)=((rhoji1(j,i,1))^2+(rhoji1(j,i,2))^2)-(((lhat(i,1)*xji1(j,i))+(lhat(i,2)*yji1(j,i)))^2);
        %--7 R1
        R1(j,i)= sqrt(((rhoji1(j,i,1))^2+(rhoji1(j,i,2))^2));
        %--8
        RN(j,i)= ((xji1(j,i))*(nhat(i,1)))+((yji1(j,i))*(nhat(i,2)));
        %--9
        RL(j,i)=(lhat(i,1)*xji1(j,i))+(lhat(i,2)*yji1(j,i));
        %--10
        L1(j,i)= ((lhat(i,1))*(nhat(i,1)))+((lhat(i,2))*(nhat(i,2)));
        %--11
        R18(j,i)=((rhoji1(j,i,1))^2+(rhoji1(j,i,2))^2)+(R12(j,i))^2-2*(RL(j,i))*(R12(j,i));
    end 
end

for j=Ns+1:N
   for i=1:N
      if j==i
          I(j,i)=(0.5*(L1(j,i))*log(((R1(j,i))^2)/R18(j,i)))+(L1(j,i))*(RL(j,i))-(RN(j,i))*((1/RL(j,i))+(1/(R12(j,i)-RL(j,i))));
      else
          I(j,i)=(0.5*(L1(j,i))*log(((R1(j,i))^2)/R18(j,i)))+(L1(j,i))*(RL(j,i))-(RN(j,i))*((1/sqrt(delta(j,i)))*(atan((RL(j,i))*((1/sqrt(delta(j,i)))))-atan(((R12(j,i))-(RL(j,i)))*((1/sqrt(delta(j,i)))))));
      end
   end
end

%% ------------------------------ Deriving V Matrix -----------------------
%% ---------------------------------- Eq.(7) ------------------------------
%%   j=1:cd1 --> x=13-13.5 , y=1.1
%    cd1
for j=1:cd1
    V(j,1)=Vp;
    for n=1:N
        Z(j,n)=(1/2*pi*eps0)*(deltat*log(K) + (Rj2(j,n)-Rj1(j,n) + (2*pi*eps0)*(R12(j,n) - R18(j,n) * atan(R18(j,n)*R12(j,n)/S(j,n))))); 
    end
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+1):(cd1+cd2) --> x=13.75-16.75 , y=0.1
%    cd2
for j=(cd1+1):(cd1+cd2)
    V(j,1)=Vn;
    for n=1:N
        Z(j,n)=(1/2*pi*eps0)*(deltat*log(K) + (Rj2(j,n)-Rj1(j,n) + (2*pi*eps0)*(R12(j,n) - R18(j,n) * atan(R18(j,n)*R12(j,n)/S(j,n))))); 
    end
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+1):(cd1+cd2+cd3) --> x=16.75 , y=0:0.1
%    cd3
for j=(cd1+cd2+1):(cd1+cd2+cd3)
    V(j,1)=Vn;
    for n=1:N
        Z(j,n)=(1/2*pi*eps0)*(deltat*log(K) + (Rj2(j,n)-Rj1(j,n) + (2*pi*eps0)*(R12(j,n) - R18(j,n) * atan(R18(j,n)*R12(j,n)/S(j,n))))); 
    end
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+cd3+1):(cd1+cd2+cd3+cd4) --> x=13.75 , y=0:0.1
%    cd4
for j=(cd1+cd2+cd3+1):(cd1+cd2+cd3+cd4)
    V(j,1)=Vn;
    for n=1:N
        Z(j,n)=(1/2*pi*eps0)*(deltat*log(K) + (Rj2(j,n)-Rj1(j,n) + (2*pi*eps0)*(R12(j,n) - R18(j,n) * atan(R18(j,n)*R12(j,n)/S(j,n))))); 
    end
end
%    Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+cd3+cd4+1):(cd1+cd2+cd3+cd4+cd5) --> x=13.5 , y=1.1:1.2
%    cd5
for j=(cd1+cd2+cd3+cd4+1):(cd1+cd2+cd3+cd4+cd5)
    V(j,1)=Vp;
    for n=1:N
        Z(j,n)=(1/2*pi*eps0)*(deltat*log(K) + (Rj2(j,n)-Rj1(j,n) + (2*pi*eps0)*(R12(j,n) - R18(j,n) * atan(R18(j,n)*R12(j,n)/S(j,n))))); 
    end
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+cd3+cd4ca1+1):(cd1+cd2+cd3+cd4+cd5+ca1) --> x=13:13.5 , y=1.2
%    ca1
for j=(cd1+cd2+cd3+cd4+ca1+1):(cd1+cd2+cd3+cd4+cd5+ca1)
    V(j,1)=Vp;
    for n=1:N
        Z(j,n)=(1/2*pi*eps0)*(deltat*log(K) + (Rj2(j,n)-Rj1(j,n) + (2*pi*eps0)*(R12(j,n) - R18(j,n) * atan(R18(j,n)*R12(j,n)/S(j,n))))); 
    end
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+cd3+cd4ca1+ca2+1):Ns --> x=13 , y=1.1:1.2
%    ca2
for j=(cd1+cd2+cd3+cd4+ca1+ca2+1):Ns
    V(j,1)=Vp;
    for n=1:N
        Z(j,n)=(1/2*pi*eps0)*(deltat*log(K) + (Rj2(j,n)-Rj1(j,n) + (2*pi*eps0)*(R12(j,n) - R18(j,n) * atan(R18(j,n)*R12(j,n)/S(j,n))))); 
    end
end


%% ---------------------------- Eq.(18) -----------------------------------
%%   j=(Ns+1):Ns+da1
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da1
for j=(Ns+1):Ns+da1
    V(j,1)=0;
    for n=1:N
        if j==n
        Z(j,j)=((eps0+eps2)/(eps0))+(-((eps2-eps0)/(2*pi*eps0)))*I(j,j);
        else
        Z(j,n)=(-((eps2-eps0)/(2*pi*eps0)))*I(j,n);
        end
    end
end

%%   j=(Ns+1)+da1:Ns+da1+da2
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da2
for j=(Ns+1)+da1:Ns+da1+da2
    V(j,1)=0;
    for n=1:N
        if j==n
        Z(j,j)=((eps0+eps2)/(eps0))+(-((eps2-eps0)/(2*pi*eps0)))*I(j,j);
        else
        Z(j,n)=(-((eps2-eps0)/(2*pi*eps0)))*I(j,n);
        end
    end
end
%%   j=(Ns+1)+da1+da2:Ns+da1+da2+da3
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da3
for j=(Ns+1)+da1+da2:Ns+da1+da2+da3
    V(j,1)=0;
    for n=1:N
        if j==n
        Z(j,j)=((eps0+eps3)/(eps0))+(-((eps3-eps0)/(2*pi*eps0)))*I(j,j);
        else
        Z(j,n)=(-((eps3-eps0)/(2*pi*eps0)))*I(j,n);
        end
    end
end
%%   j=(Ns+1)+da1+da2+da3:Ns+da1+da2+da3+da4
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da4
for j=(Ns+1)+da1+da2+da3:Ns+da1+da2+da3+da4
    V(j,1)=0;
    for n=1:N
        if j==n
        Z(j,j)=((eps0+eps3)/(eps0))+(-((eps3-eps0)/(2*pi*eps0)))*I(j,j);
        else
        Z(j,n)=(-((eps3-eps0)/(2*pi*eps0)))*I(j,n);
        end
    end
end
%%   j=(Ns+1)+da1++da2+da3+da4:Ns+da1+da2+da3+da4+da5
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da5
for j=(Ns+1)+da1++da2+da3+da4:Ns+da1+da2+da3+da4+da5
    V(j,1)=0;
    for n=1:N
        if j==n
        Z(j,j)=((eps0+eps1)/(eps0))+(-((eps1-eps0)/(2*pi*eps0)))*I(j,j);
        else
        Z(j,n)=(-((eps1-eps0)/(2*pi*eps0)))*I(j,n);
        end
    end
end
%%   j=(Ns+1)+da1+da2+da3+da4+da5:Ns+da1+da2+da3+da4+da5+dd1
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd1
for j=(Ns+1)+da1+da2+da3+da4+da5:Ns+da1+da2+da3+da4+da5+dd1
    V(j,1)=0;
    for n=1:N
        if j==n
        Z(j,j)=((eps3+eps1)/(eps0))+(-((eps1-eps3)/(2*pi*eps0)))*I(j,j);
        else
        Z(j,n)=(-((eps1-eps3)/(2*pi*eps0)))*I(j,n);
        end
    end
end
%%   j=(Ns+1)+da1+da2+da3+da4+da5+dd1:Ns+da1+da2+da3+da4+da5+dd1+dd2
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd2
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1:Ns+da1+da2+da3+da4+da5+dd1+dd2
    V(j,1)=0;
    for n=1:N
        if j==n
        Z(j,j)=((eps1+eps2)/(eps0))+(-((eps1-eps2)/(2*pi*eps0)))*I(j,j);
        else
        Z(j,n)=(-((eps1-eps2)/(2*pi*eps0)))*I(j,n);
        end
    end
end
%%  j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2:Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd3
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2:Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
    V(j,1)=0;
    for n=1:N
        if j==n
        Z(j,j)=((eps1+eps2)/(eps0))+(-((eps1-eps2)/(2*pi*eps0)))*I(j,j);
        else
        Z(j,n)=(-((eps1-eps2)/(2*pi*eps0)))*I(j,n);
        end
    end
end
%%   j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2+dd3:N
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd4
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2+dd3:N
    V(j,1)=0;
    for n=1:N
        if j==n
        Z(j,j)=((eps1+eps2)/(eps0))+(-((eps1-eps2)/(2*pi*eps0)))*I(j,j);
        else
        Z(j,n)=(-((eps1-eps2)/(2*pi*eps0)))*I(j,n);
        end
    end
end
%% ------------------- Trying To Solve ------------------------------------
sigma(:,1)= linsolve(Z,V);


