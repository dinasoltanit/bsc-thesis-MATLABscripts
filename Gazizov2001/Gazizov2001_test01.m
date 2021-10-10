clc;
clear;
%syms taw
% this is not symbolic
%-----------------------------General-------------------------------------**
ref=0.025; Vn=0; Vp=1; %Vm=4000*sqrt(2); Vp=Vm*sin(2*pi*fm*taw);
eps0=8.85*10^(-12); eps1=1.5;   eps2=4; eps3=(1-25.935j);  %eps of plasma which is a function
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
%gnd=l6/ref;
%-------------------------------- ground ----------------------------------
Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2; %without gnd-to-conductor interface
Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
N= Ns + Nd;
%-------------------------------- Coordinate Matrix ----------------------**
Z=zeros(N,N);
V=(zeros(N,1));
cdensity=(zeros(N,1));
temp=(zeros(N,1));
Q=(zeros(16,1));
Density_q=(zeros(16,1));
%% ------------------------------ Dami Matrices ---------------------------
r=zeros(N,5); %rs=[j,xj,yj,dj,thetaj]
a1=zeros(N,N);
b1=zeros(N,N);
a2=zeros(N,N);
b2=zeros(N,N);
I=zeros(N,N);
Ihat=zeros(N,N);
%% ----------------------------- Deriving rs Matrix ----------------------**
tic;
% cd1
for j=1:cd1
    %rs=[j,xj,yj]
    r(j,1)=j;
    r(j,2)=13+(ref/2)+(j-1)*ref; r(j,3)=1.1; r(j,4)=ref; r(j,5)=0;
end
% cd2
for j=1+cd1:cd2+cd1
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13.75+(ref/2)+(j-1)*ref; r(j,3)=0.1; r(j,4)=ref; r(j,5)=pi;
end
%cd3
for j=1+cd1+cd2:cd3+cd1+cd2
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=16.75; r(j,3)=0+(ref/2)+(j-1)*ref; r(j,4)=ref; r(j,5)=pi/2;
end
%cd4
for j=1+cd1+cd2+cd3:cd4+cd1+cd2+cd3
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13.75; r(j,3)=0+(ref/2)+(j-1)*ref; r(j,4)=ref; r(j,5)=(3*pi)/2;
end
%cd5
for j=1+cd1+cd2+cd3+cd4:cd5+cd1+cd2+cd3+cd4
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13.5; r(j,3)=1.1+(ref/2)+(j-1)*ref; r(j,4)=ref; r(j,5)=pi/2;
end
%ca1
for j=1+cd1+cd2+cd3+cd4+cd5:ca1+cd1+cd2+cd3+cd4+cd5
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13+(ref/2)+(j-1)*ref; r(j,3)=1.2; r(j,4)=ref; r(j,5)=pi;
end
%ca2
for j=1+cd1+cd2+cd3+cd4+cd5+ca1:ca2+cd1+cd2+cd3+cd4+cd5+ca1
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13; r(j,3)=1.1+(ref/2)+(j-1)*ref; r(j,4)=ref; r(j,5)=(3*pi)/2;
end
fprintf('rs matrix Found.\n')
toc;
%% ----------------------------- Deriving rd Matrix ----------------------**
tic;
%r=[i,xi,yi]
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da1
for j=1+Ns:da1+Ns
   r(j,1)=j; r(j,2)=0+(ref/2)+(j-1)*ref; r(j,3)=1.1; r(j,4)=ref; r(j,5)=pi;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da2
for j=1+Ns+da1:da2+Ns+da1
   r(j,1)=j; r(j,2)=16.75+(ref/2)+(j-1)*ref; r(j,3)=1.1; r(j,4)=ref; r(j,5)=pi;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da3
for j=1+Ns+da1+da2:da3+Ns+da1+da2
   r(j,1)=j; r(j,2)=13.5; r(j,3)=1.2+(ref/2)+(j-1)*ref; r(j,4)=ref; r(j,5)=(3*pi)/2;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da4
for j=1+Ns+da1+da2+da3:da4+Ns+da1+da2+da3
    r(j,1)=j; r(j,2)=13.5+(b/(2*da4))+(j-1)*(b/da4); r(j,3)=1.2+(a/(2*da4))+(j-1)*(a/da4); r(j,4)=(sqrt(a^2+b^2))/(da4); r(j,5)=pi-atan(a/b);
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da5
for j=1+Ns+da1+da2+da3+da4:da5+Ns+da1+da2+da3+da4
   r(j,1)=j; r(j,2)=16.5+(ref/2)+(j-1)*ref; r(j,3)=1.1; r(j,4)=ref; r(j,5)=pi;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd1
for j=1+Ns+da1+da2+da3+da4+da5:dd1+Ns+da1+da2+da3+da4+da5
   r(j,1)=j; r(j,2)=13.5+(ref/2)+(j-1)*ref; r(j,3)=1.1; r(j,4)=ref; r(j,5)=pi;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd2
for j=1+Ns+da1+da2+da3+da4+da5+dd1:dd2+Ns+da1+da2+da3+da4+da5+dd1
   r(j,1)=j; r(j,2)=16.75; r(j,3)=0.1+(ref/2)+(j-1)*ref; r(j,4)=ref; r(j,5)=pi/2;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd3
for j=1+Ns+da1+da2+da3+da4+da5+dd1+dd2:dd3+Ns+da1+da2+da3+da4+da5+dd1+dd2
   r(j,1)=j; r(j,2)=13; r(j,3)=0.1+(ref/2)+(j-1)*ref; r(j,4)=ref; r(j,5)=(3*pi)/2;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd4
for j=1+Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3:dd4+Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
   r(j,1)=j; r(j,2)=13+(ref/2)+(j-1)*ref; r(j,3)=0.1; r(j,4)=ref; r(j,5)=pi;
end
toc
fprintf('rd Matrix Found.\n')
%% -------------------- Finding a1 b1 a2 b2 I Ihat ------------------------
for m=1:N
   for n=1:N
      % r=zeros(N,5); r=[j,xj,yj,dj,thetaj]
        a1(m,n)=((r(m,2))-(r(n,2)))*sin(r(n,5)) + ((r(m,3))+(r(n,3)))*cos(r(n,5));
        b1(m,n)=((r(m,2))-(r(n,2)))*cos(r(n,5)) - ((r(m,3))+(r(n,3)))*sin(r(n,5));
        a2(m,n)=((r(m,2))-(r(n,2)))*sin(r(n,5)) - ((r(m,3))-(r(n,3)))*cos(r(n,5)); 
        b2(m,n)=((r(m,2))-(r(n,2)))*cos(r(n,5)) + ((r(m,3))-(r(n,3)))*sin(r(n,5));
   end   
end

%% -------------------------- Finding Z Matrix ----------------------------
tic;
% for 1:Ns
for m=1:Ns
   for n=1:N
        Z(m,n)=(1/4*pi*eps0)*(F1(a1(m,n),b1(m,n),r(n,4)) - F1(a2(m,n),b2(m,n),r(n,4)));
   end    
end
% for da1 da2 da5 dd1 dd4
%da1 , da2
for m=1+Ns:da1+da2+Ns
   for n=1:N
       % r=zeros(N,5); r=[j,xj,yj,dj,thetaj]
       I(m,n)=(r(m,3)-r(n,3)-(b2(m,n))*sin(r(n,5)))*F2(a2(m,n),b2(m,n),r(n,4))-(sin(r(n,5)))*F3(a2(m,n),b2(m,n),r(n,4));
       Ihat(m,n)=(r(m,3)+r(n,3)+(b2(m,n))*sin(r(n,5)))*F2(a1(m,n),b1(m,n),r(n,4))-(sin(r(n,5)))*F3(a1(m,n),b1(m,n),r(n,4));
       if m==n
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((1+eps2)/(1-eps2));
       else
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n));
       end
   end
end
%da5
for m=1+Ns+da1+da2+da3+da4:Ns+da1+da2+da3+da4+da5
   for n=1:N
       % r=zeros(N,5); r=[j,xj,yj,dj,thetaj]
       I(m,n)=(r(m,3)-r(n,3)-(b2(m,n))*sin(r(n,5)))*F2(a2(m,n),b2(m,n),r(n,4))-(sin(r(n,5)))*F3(a2(m,n),b2(m,n),r(n,4));
       Ihat(m,n)=(r(m,3)+r(n,3)+(b2(m,n))*sin(r(n,5)))*F2(a1(m,n),b1(m,n),r(n,4))-(sin(r(n,5)))*F3(a1(m,n),b1(m,n),r(n,4));
       if m==n
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((1+eps1)/(1-eps1));
       else
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n));
       end
   end   
end
%for dd1
for m=1+Ns+da1+da2+da3+da4+da5:Ns+da1+da2+da3+da4+da5+dd1
   for n=1:N
       % r=zeros(N,5); r=[j,xj,yj,dj,thetaj]
       I(m,n)=(r(m,3)-r(n,3)-(b2(m,n))*sin(r(n,5)))*F2(a2(m,n),b2(m,n),r(n,4))-(sin(r(n,5)))*F3(a2(m,n),b2(m,n),r(n,4));
       Ihat(m,n)=(r(m,3)+r(n,3)+(b2(m,n))*sin(r(n,5)))*F2(a1(m,n),b1(m,n),r(n,4))-(sin(r(n,5)))*F3(a1(m,n),b1(m,n),r(n,4));
       if m==n
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((eps3+eps1)/(eps3-eps1));
       else
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n));
       end
   end   
end
%for dd4
for m=1+Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3:N
   for n=1:N
       % r=zeros(N,5); r=[j,xj,yj,dj,thetaj]
       I(m,n)=(r(m,3)-r(n,3)-(b2(m,n))*sin(r(n,5)))*F2(a2(m,n),b2(m,n),r(n,4))-(sin(r(n,5)))*F3(a2(m,n),b2(m,n),r(n,4));
       Ihat(m,n)=(r(m,3)+r(n,3)+(b2(m,n))*sin(r(n,5)))*F2(a1(m,n),b1(m,n),r(n,4))-(sin(r(n,5)))*F3(a1(m,n),b1(m,n),r(n,4));
       if m==n
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((eps1+eps2)/(eps1-eps2));
       else
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n));
       end
   end   
end
%for da3 dd2 dd3
%dd2 , dd3
for m=1+Ns+da1+da2+da3+da4+da5+dd1:Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
   for n=1:N
       % r=zeros(N,5); r=[j,xj,yj,dj,thetaj]
       I(m,n)=(r(m,2)-r(n,2)-(b2(m,n))*cos(r(n,5)))*F2(a2(m,n),b2(m,n),r(n,4))-(cos(r(n,5)))*F3(a2(m,n),b2(m,n),r(n,4));
       Ihat(m,n)=(r(m,2)-r(n,2)-(b2(m,n))*cos(r(n,5)))*F2(a1(m,n),b1(m,n),r(n,4))-(cos(r(n,5)))*F3(a1(m,n),b1(m,n),r(n,4));
       if m==n
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((eps2+eps1)/(eps2-eps1));
       else
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n));
       end
   end
end
% for da3
for m=1+Ns+da1+da2:da3+Ns+da1+da2
   for n=1:N
       % r=zeros(N,5); r=[j,xj,yj,dj,thetaj]
       I(m,n)=(r(m,2)-r(n,2)-(b2(m,n))*cos(r(n,5)))*F2(a2(m,n),b2(m,n),r(n,4))-(cos(r(n,5)))*F3(a2(m,n),b2(m,n),r(n,4));
       Ihat(m,n)=(r(m,2)-r(n,2)-(b2(m,n))*cos(r(n,5)))*F2(a1(m,n),b1(m,n),r(n,4))-(cos(r(n,5)))*F3(a1(m,n),b1(m,n),r(n,4));
       if m==n
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((1+eps3)/(1-eps3));
       else
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n));
       end
   end   
end
% for da4
for m=1+Ns+da1+da2+da3:da4+Ns+da1+da2+da3
   for n=1:N
       % r=zeros(N,5); r=[j,xj,yj,dj,thetaj]
       I(m,n)=((sin(r(m,4)))*((r(m,2)-r(n,2)-(b2(m,n))*cos(r(n,5)))*F2(a2(m,n),b2(m,n),r(n,4))-(cos(r(n,5)))*F3(a2(m,n),b2(m,n),r(n,4))))-((cos(r(m,4)))*((r(m,3)-r(n,3)-(b2(m,n))*sin(r(n,5)))*F2(a2(m,n),b2(m,n),r(n,4))-(sin(r(n,5)))*F3(a2(m,n),b2(m,n),r(n,4))));
    Ihat(m,n)=((sin(r(m,4))))*((r(m,2)-r(n,2)-(b2(m,n))*cos(r(n,5)))*F2(a1(m,n),b1(m,n),r(n,4))-(cos(r(n,5)))*F3(a1(m,n),b1(m,n),r(n,4)))-(cos(r(m,4)))*((r(m,3)+r(n,3)+(b2(m,n))*sin(r(n,5)))*F2(a1(m,n),b1(m,n),r(n,4))-(sin(r(n,5)))*F3(a1(m,n),b1(m,n),r(n,4)));
       if m==n
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((1+eps3)/(1-eps3));
       else
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n));
       end
   end    
end
toc
fprintf('Z matrix Found.\n')
%% -------------------------- Finding V Matrix ----------------------------
tic;
%%   j=1:cd1 --> x=13-13.5 , y=1.1
%    cd1
for j=1:cd1
    V(j,1)=Vp;
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+1):(cd1+cd2) --> x=13.75-16.75 , y=0.1
%    cd2
for j=(cd1+1):(cd1+cd2)
    V(j,1)=Vn;
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+1):(cd1+cd2+cd3) --> x=16.75 , y=0:0.1
%    cd3
for j=(cd1+cd2+1):(cd1+cd2+cd3)
    V(j,1)=Vn;
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+cd3+1):(cd1+cd2+cd3+cd4) --> x=13.75 , y=0:0.1
%    cd4
for j=(cd1+cd2+cd3+1):(cd1+cd2+cd3+cd4)
    V(j,1)=Vn;
end
%    Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+cd3+cd4+1):(cd1+cd2+cd3+cd4+cd5) --> x=13.5 , y=1.1:1.2
%    cd5
for j=(cd1+cd2+cd3+cd4+1):(cd1+cd2+cd3+cd4+cd5)
    V(j,1)=Vp;
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+cd3+cd4ca1+1):(cd1+cd2+cd3+cd4+cd5+ca1) --> x=13:13.5 , y=1.2
%    ca1
for j=(cd1+cd2+cd3+cd4+ca1+1):(cd1+cd2+cd3+cd4+cd5+ca1)
    V(j,1)=Vp;
end
% Ns= cd1+cd2+cd3+cd4+cd5+ca1+ca2;
%%   j=(cd1+cd2+cd3+cd4ca1+ca2+1):Ns --> x=13 , y=1.1:1.2
%    ca2
for j=(cd1+cd2+cd3+cd4+ca1+ca2+1):Ns
    V(j,1)=Vp;
end


%% ---------------------------- Eq.(18) -----------------------------------
%%   j=(Ns+1):Ns+da1
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da1
for j=(Ns+1):Ns+da1
    V(j,1)=0;
end

%%   j=(Ns+1)+da1:Ns+da1+da2
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da2
for j=(Ns+1)+da1:Ns+da1+da2
    V(j,1)=0;
end
%%   j=(Ns+1)+da1+da2:Ns+da1+da2+da3
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da3
for j=(Ns+1)+da1+da2:Ns+da1+da2+da3
    V(j,1)=0;
end
%%   j=(Ns+1)+da1+da2+da3:Ns+da1+da2+da3+da4
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da4
for j=(Ns+1)+da1+da2+da3:Ns+da1+da2+da3+da4
    V(j,1)=0;
end
%%   j=(Ns+1)+da1++da2+da3+da4:Ns+da1+da2+da3+da4+da5
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da5
for j=(Ns+1)+da1++da2+da3+da4:Ns+da1+da2+da3+da4+da5
    V(j,1)=0;
end
%%   j=(Ns+1)+da1+da2+da3+da4+da5:Ns+da1+da2+da3+da4+da5+dd1
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd1
for j=(Ns+1)+da1+da2+da3+da4+da5:Ns+da1+da2+da3+da4+da5+dd1
    V(j,1)=0;
end
%%   j=(Ns+1)+da1+da2+da3+da4+da5+dd1:Ns+da1+da2+da3+da4+da5+dd1+dd2
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd2
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1:Ns+da1+da2+da3+da4+da5+dd1+dd2
    V(j,1)=0;
end
%%  j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2:Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd3
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2:Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
    V(j,1)=0;
end
%%   j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2+dd3:N
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd4
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2+dd3:N
    V(j,1)=0;
end
toc
fprintf('V matrix Found.\n')
%% ---------------------- Finding Charge Densities ------------------------
cond(Z)
cdensity(:,1)= linsolve(Z,V);
%% -------------------------- Finding Charges -----------------------------
%    cd1
for j=1:cd1
    temp(j,1)=0;
    Q(1,1)= temp(j,1);
    Q(1,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(1,1)=(Q(1,1))/l1;
end
%    cd2
for j=(cd1+1):(cd1+cd2)
    temp(j,1)=0;
    Q(2,1)= temp(j,1);
    Q(2,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(2,1)=(Q(2,1))/l2;
end
%    cd3
for j=(cd1+cd2+1):(cd1+cd2+cd3)
    temp(j,1)=0;
    Q(3,1)= temp(j,1);
    Q(3,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(3,1)=(Q(3,1))/h2;
end
%    cd4
for j=(cd1+cd2+cd3+1):(cd1+cd2+cd3+cd4)
    temp(j,1)=0;
    Q(4,1)= temp(j,1);
    Q(4,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(4,1)=(Q(4,1))/h2;
end
%    cd5
for j=(cd1+cd2+cd3+cd4+1):(cd1+cd2+cd3+cd4+cd5)
    temp(j,1)=0;
    Q(5,1)= temp(j,1);
    Q(5,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(5,1)=(Q(5,1))/h1;
end
%    ca1
for j=(cd1+cd2+cd3+cd4+ca1+1):(cd1+cd2+cd3+cd4+cd5+ca1)
    temp(j,1)=0;
    Q(6,1)= temp(j,1);
    Q(6,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(6,1)=(Q(6,1))/l1;
end
%    ca2
for j=(cd1+cd2+cd3+cd4+ca1+ca2+1):Ns
    temp(j,1)=0;
    Q(7,1)= temp(j,1);
    Q(7,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(7,1)=(Q(7,1))/h1;
end
%da1
for j=(Ns+1):Ns+da1
    temp(j,1)=0;
    Q(8,1)= temp(j,1);
    Q(8,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(8,1)=(Q(8,1))/l8;
end
%da2
for j=(Ns+1)+da1:Ns+da1+da2
    temp(j,1)=0;
    Q(9,1)= temp(j,1);
    Q(9,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(9,1)=(Q(9,1))/l9;
end
%da3
for j=(Ns+1)+da1+da2:Ns+da1+da2+da3
    temp(j,1)=0;
    Q(10,1)= temp(j,1);
    Q(10,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(10,1)=(Q(10,1))/a;
end
%da4
for j=(Ns+1)+da1+da2+da3:Ns+da1+da2+da3+da4
    temp(j,1)=0;
    Q(11,1)= temp(j,1);
    Q(11,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(11,1)=(Q(11,1))/a;
end
%da5
for j=(Ns+1)+da1++da2+da3+da4:Ns+da1+da2+da3+da4+da5
    temp(j,1)=0;
    Q(12,1)= temp(j,1);
    Q(12,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(12,1)=(Q(12,1))/l4;
end
%dd1
for j=(Ns+1)+da1+da2+da3+da4+da5:Ns+da1+da2+da3+da4+da5+dd1
    temp(j,1)=0;
    Q(13,1)= temp(j,1);
    Q(13,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(13,1)=(Q(13,1))/b;
end
%dd2
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1:Ns+da1+da2+da3+da4+da5+dd1+dd2
    temp(j,1)=0;
    Q(14,1)= temp(j,1);
    Q(14,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(14,1)=(Q(14,1))/h4;
end
%dd3
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2:Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
    temp(j,1)=0;
    Q(15,1)= temp(j,1);
    Q(15,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(15,1)=(Q(15,1))/h4;
end
%dd4
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2+dd3:N
    temp(j,1)=0;
    Q(16,1)= temp(j,1);
    Q(16,1)= temp(j,1) + (cdensity(j,1))*r(j,5);
    Density_q(16,1)=(Q(16,1))/l5;
end

fprintf('All Done! \n')