clc;
clear;
%syms taw
% this is not symbolic
%-----------------------------General-------------------------------------**
ref=0.025e-3; Vn=0; Vp=1; %Vm=4000*sqrt(2); Vp=Vm*sin(2*pi*fm*taw);
eps0=8.85*10^(-12); eps1=2;   %eps2=4; eps3=(1-25.935j);  %eps of plasma which is a function
%-----------------------------Dimensions----------------------------------*
h1=1e-3; %mm
l1=2e-3;
%--
h2=1e-3; %mm
l2=3e-3;
%--
h3=1e-3; %mm
l3=3e-3;
%--
h4=1e-3; %mm
l4=3e-3;
%--
h5=1e-3; %mm
l5=3e-3;
%--
h6=1e-3; %mm
l6=6e-3;
%--
h7=1e-3; %mm
l7=6e-3;
%--------------------- dielectric-to-dielectric interfaces ----------------
da1=10; da2=8 ; da3=10;
%--------------------- dielectric-to-conductor interfaces -----------------
ca1=2; ca3=2; ca4=2; ca6=2;
ca2=8; ca5=8;
cd1=8; cd2=8;
%-------------------------------- ground ----------------------------------
%gnd=l6/ref;
%-------------------------------- ground ----------------------------------
Ns= cd1+cd2+ca1+ca2+ca3+ca4+ca5+ca6; %without gnd-to-conductor interface
Nd= da1+da2+da3;
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
    r(j,2)=6*(10^(-3))+(ref/2)+(j-1)*ref; r(j,3)=1*(10^(-3)); r(j,4)=l3/cd1; r(j,5)=0;
end
% cd2
for j=1+cd1:cd2+cd1
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=11*(10^(-3))+(ref/2)+(j-1)*ref; r(j,3)=1*(10^(-3)); r(j,4)=l5/cd2; r(j,5)=0;
end
%ca1
for j=1+cd1+cd2:ca1+cd1+cd2
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=6*(10^(-3)); r(j,3)=1*(10^(-3))+(ref/2)+(j-1)*ref; r(j,4)=h3/ca1; r(j,5)=(3*pi)/2;
end
%ca2
for j=1+cd1+cd2+ca1:ca2+cd1+cd2+ca1
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=6*(10^(-3))+(ref/2)+(j-1)*ref; r(j,3)=2*(10^(-3)); r(j,4)=l3/ca2; r(j,5)=pi;
end
%ca3
for j=1+ca2+cd1+cd2+ca1:ca2+cd1+cd2+ca1+ca3
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=9*(10^(-3)); r(j,3)=1*(10^(-3))+(ref/2)+(j-1)*ref; r(j,4)=h3/ca3; r(j,5)=pi/2;
end
%ca4
for j=1+ca2+cd1+cd2+ca1+ca3:ca2+cd1+cd2+ca1+ca3+ca4
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=11*(10^(-3)); r(j,3)=1*(10^(-3))+(ref/2)+(j-1)*ref; r(j,4)=h2/ca4; r(j,5)=(3*pi);
end
%ca5
for j=1+ca2+cd1+cd2+ca1+ca3+ca4:ca2+cd1+cd2+ca1+ca3+ca4+ca5
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=11*(10^(-3))+(ref/2)+(j-1)*ref; r(j,3)=2*(10^(-3)); r(j,4)=l2/ca5; r(j,5)=pi;
end
%ca6
for j=1+ca2+cd1+cd2+ca1+ca3+ca4+ca5:ca2+cd1+cd2+ca1+ca3+ca4+ca5+ca6
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=14*(10^(-3)); r(j,3)=1*(10^(-3))+(ref/2)+(j-1)*ref; r(j,4)=h2/ca6; r(j,5)=pi/2;
end
fprintf('rs matrix Found.\n')
toc;
%% ----------------------------- Deriving rd Matrix ----------------------**
tic;
%r=[i,xi,yi]
%da1
for j=1+Ns:da1+Ns
   r(j,1)=j; r(j,2)=0+(ref/2)+(j-1)*ref; r(j,3)=1*(10^(-3)); r(j,4)=l6/da1; r(j,5)=pi;
end
%da2
for j=1+Ns+da1:da2+Ns+da1
   r(j,1)=j; r(j,2)=9*(10^(-3))+(ref/2)+(j-1)*ref; r(j,3)=1*(10^(-3)); r(j,4)=l1/da2; r(j,5)=pi;
end
%da3
for j=1+Ns+da1+da2:da3+Ns+da1+da2
   r(j,1)=j; r(j,2)=14*(10^(-3))+(ref/2)+(j-1)*ref; r(j,3)=1*(10^(-3)); r(j,4)=l7/da3; r(j,5)=pi;
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
% for da1 da2 da3
%da1 , da2
for m=1+Ns:da1+da2+da3+Ns
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
%da5
%for dd1
%for dd4
%for da3 dd2 dd3
%dd2 , dd3
% for da3
% for da4
toc
fprintf('Z matrix Found.\n')
%% -------------------------- Finding V Matrix ----------------------------
tic;
%    cd1
for j=1:cd1
    V(j,1)=Vp;
end
%    cd2
for j=(cd1+1):(cd1+cd2)
    V(j,1)=Vp;
end
%    ca1
for j=(cd1+cd2+1):(cd1+cd2+ca1)
    V(j,1)=Vp;
end
%    ca2
for j=(cd1+cd2+ca1+1):cd1+cd2+ca1+ca2
    V(j,1)=Vp;
end
%    ca3
for j=(cd1+cd2+ca1+ca2+1):cd1+cd2+ca1+ca2+ca3
    V(j,1)=Vp;
end
%    ca4
for j=(cd1+cd2+ca1+ca2+ca3+1):cd1+cd2+ca1+ca2+ca3+ca4
    V(j,1)=Vp;
end
%    ca5
for j=(cd1+cd2+ca1+ca2+ca3+ca4+1):cd1+cd2+ca1+ca2+ca3+ca4+ca5
    V(j,1)=Vp;
end
%    ca6
for j=(cd1+cd2+ca1+ca2+ca3+ca4+ca5+1):Ns
    V(j,1)=Vp;
end

%% ---------------------------- Eq.(18) -----------------------------------
%da1
for j=(Ns+1):Ns+da1
    V(j,1)=0;
end
%da2
for j=(Ns+1)+da1:Ns+da1+da2
    V(j,1)=0;
end
%da3
for j=(Ns+1)+da1+da2:Ns+da1+da2+da3
    V(j,1)=0;
end
toc
fprintf('V matrix Found.\n')
%% ---------------------- Finding Charge Densities ------------------------
cond(Z)
cdensity(:,1)= linsolve(Z,V);

fprintf('All Done! \n')