clc;
clear;
%syms taw
% this is not symbolic
% in this script, I tried to evaluate charge distribution on each section.
% I don't know if I want to calculate the average charge density as it
% seems not very accurate. However, I don't know how should I give the
% calculated charge density on local section to the cfd tool.
%-----------------------------General-------------------------------------**
fm=3000;
taw=linspace(0,0.000670,100);
Nt=numel(taw);
%-----------------------------General-------------------------------------**
ref=0.05e-3; nh=(11.1e-3)/ref; nl=(21.5e-3)/ref; dx=ref; dy=ref; Vn=0; Vm=4000*sqrt(2); Vp=Vm*sin(2*pi*fm*taw);
eps0=8.85*10^(-12); eps1=1.5*eps0;   eps2=4*eps0; eps3=(1-25.935j)*eps0;  %eps of plasma which is a function
deltat=ref; K=1;
%-----------------------------Dimensions----------------------------------*
h1=0.1e-3; %mm
l1=0.5e-3;
%--
h2=0.1e-3; %mm
l2=3e-3;
%--
a=1.5e-3;
h3=a-h1;
b=3e-3; %mm %h3,l3
%--
h4=1e-3; %mm
l4=0.25e-3;
%--
h5=0.1e-3; %mm
l5=0.75e-3;
%--
l6=21.5e-3;
%--
h8=1.1e-3; %mm
l8=13e-3;
%--
h9=1.1e-3; %mm
l9=4.75e-3;

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
Z=zeros(N,N);
V=(zeros(N,Nt));
cdensity=(zeros(N,Nt));
temp=(zeros(N,Nt));
Q=(zeros(16,Nt));
Density_q=(zeros(16,Nt));
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
    r(j,2)=13+(ref/2)+(j-1)*0.05; r(j,3)=1.1; r(j,4)=ref; r(j,5)=0;
end
% cd2
for j=1+cd1:cd2+cd1
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13.75+(ref/2)+(j-1)*0.05; r(j,3)=0.1; r(j,4)=ref; r(j,5)=pi;
end
%cd3
for j=1+cd1+cd2:cd3+cd1+cd2
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=16.75; r(j,3)=0+(ref/2)+(j-1)*0.05; r(j,4)=ref; r(j,5)=pi/2;
end
%cd4
for j=1+cd1+cd2+cd3:cd4+cd1+cd2+cd3
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13.75; r(j,3)=0+(ref/2)+(j-1)*0.05; r(j,4)=ref; r(j,5)=(3*pi)/2;
end
%cd5
for j=1+cd1+cd2+cd3+cd4:cd5+cd1+cd2+cd3+cd4
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13.5; r(j,3)=1.1+(ref/2)+(j-1)*0.05; r(j,4)=ref; r(j,5)=pi/2;
end
%ca1
for j=1+cd1+cd2+cd3+cd4+cd5:ca1+cd1+cd2+cd3+cd4+cd5
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13+(ref/2)+(j-1)*0.05; r(j,3)=1.2; r(j,4)=ref; r(j,5)=pi;
end
%ca2
for j=1+cd1+cd2+cd3+cd4+cd5+ca1:ca2+cd1+cd2+cd3+cd4+cd5+ca1
    %r=[j,xj,yj]
    r(j,1)=j; r(j,2)=13; r(j,3)=1.1+(ref/2)+(j-1)*0.05; r(j,4)=ref; r(j,5)=(3*pi)/2;
end
fprintf('rs matrix Found.\n')
toc;
%% ----------------------------- Deriving rd Matrix ----------------------**
tic;
%r=[i,xi,yi]
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da1
for j=1+Ns:da1+Ns
   r(j,1)=j; r(j,2)=0+(ref/2)+(j-1)*0.05; r(j,3)=1.1; r(j,4)=ref; r(j,5)=pi;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da2
for j=1+Ns+da1:da2+Ns+da1
   r(j,1)=j; r(j,2)=16.75+(ref/2)+(j-1)*0.05; r(j,3)=1.1; r(j,4)=ref; r(j,5)=pi;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da3
for j=1+Ns+da1+da2:da3+Ns+da1+da2
   r(j,1)=j; r(j,2)=13.5; r(j,3)=1.2+(ref/2)+(j-1)*0.05; r(j,4)=ref; r(j,5)=(3*pi)/2;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da4
for j=1+Ns+da1+da2+da3:da4+Ns+da1+da2+da3
    r(j,1)=j; r(j,2)=13.5+(b/(2*da4))+(j-1)*(b/da4); r(j,3)=1.2+(a/(2*da4))+(j-1)*(a/da4); r(j,4)=(sqrt(a^2+b^2))/(da4); r(j,5)=pi-atan(a/b);
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%da5
for j=1+Ns+da1+da2+da3+da4:da5+Ns+da1+da2+da3+da4
   r(j,1)=j; r(j,2)=16.5+(ref/2)+(j-1)*0.05; r(j,3)=1.1; r(j,4)=ref; r(j,5)=pi;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd1
for j=1+Ns+da1+da2+da3+da4+da5:dd1+Ns+da1+da2+da3+da4+da5
   r(j,1)=j; r(j,2)=13.5+(ref/2)+(j-1)*0.05; r(j,3)=1.1; r(j,4)=ref; r(j,5)=pi;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd2
for j=1+Ns+da1+da2+da3+da4+da5+dd1:dd2+Ns+da1+da2+da3+da4+da5+dd1
   r(j,1)=j; r(j,2)=16.75; r(j,3)=0.1+(ref/2)+(j-1)*0.05; r(j,4)=ref; r(j,5)=pi/2;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd3
for j=1+Ns+da1+da2+da3+da4+da5+dd1+dd2:dd3+Ns+da1+da2+da3+da4+da5+dd1+dd2
   r(j,1)=j; r(j,2)=13; r(j,3)=0.1+(ref/2)+(j-1)*0.05; r(j,4)=ref; r(j,5)=(3*pi)/2;
end
%Nd= da1+da2+da3+da4+da5+dd1+dd2+dd3+dd4;
%dd4
for j=1+Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3:dd4+Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
   r(j,1)=j; r(j,2)=13+(ref/2)+(j-1)*0.05; r(j,3)=0.1; r(j,4)=ref; r(j,5)=pi;
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
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((eps0+eps2)/(eps0-eps2));
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
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((eps0+eps1)/(eps0-eps1));
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
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((eps0+eps3)/(eps0-eps3));
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
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n))+(1/2*eps0)*((eps0+eps3)/(eps0-eps3));
       else
           Z(m,n)=(1/2*pi*eps0)*(I(m,n)-Ihat(m,n));
       end
   end    
end
toc
fprintf('Z matrix Found.\n')
%% -------------------------- Finding V Matrix ----------------------------
for k=1:Nt
    tic;
    %    cd1
    for j=1:cd1
        V(j,k)=Vp(1,k);
    end
    %    cd2
    for j=(cd1+1):(cd1+cd2)
        V(j,k)=Vn;
    end
    %    cd3
    for j=(cd1+cd2+1):(cd1+cd2+cd3)
        V(j,k)=Vn;
    end
    %    cd4
    for j=(cd1+cd2+cd3+1):(cd1+cd2+cd3+cd4)
        V(j,k)=Vn;
    end
    %    cd5
    for j=(cd1+cd2+cd3+cd4+1):(cd1+cd2+cd3+cd4+cd5)
        V(j,k)=Vp(1,k);
    end
    %    ca1
    for j=(cd1+cd2+cd3+cd4+ca1+1):(cd1+cd2+cd3+cd4+cd5+ca1)
        V(j,k)=Vp(1,k);
    end
    %    ca2
    for j=(cd1+cd2+cd3+cd4+ca1+ca2+1):Ns
        V(j,k)=Vp(1,k);
    end

    %% ---------------------------- Eq.(18) -----------------------------------
    %da1
    for j=(Ns+1):Ns+da1
        V(j,k)=0;
    end
    %da2
    for j=(Ns+1)+da1:Ns+da1+da2
        V(j,k)=0;
    end
    %da3
    for j=(Ns+1)+da1+da2:Ns+da1+da2+da3
        V(j,k)=0;
    end
    %da4
    for j=(Ns+1)+da1+da2+da3:Ns+da1+da2+da3+da4
        V(j,k)=0;
    end
    %da5
    for j=(Ns+1)+da1++da2+da3+da4:Ns+da1+da2+da3+da4+da5
        V(j,k)=0;
    end
    %dd1
    for j=(Ns+1)+da1+da2+da3+da4+da5:Ns+da1+da2+da3+da4+da5+dd1
        V(j,k)=0;
    end
    %dd2
    for j=(Ns+1)+da1+da2+da3+da4+da5+dd1:Ns+da1+da2+da3+da4+da5+dd1+dd2
        V(j,k)=0;
    end
    %dd3
    for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2:Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
        V(j,k)=0;
    end
    %dd4
    for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2+dd3:N
        V(j,k)=0;
    end
    toc
fprintf('V matrix Found.\n')
%% ---------------------- Finding Charge Densities ------------------------
cond(Z)
cdensity(:,k)= linsolve(Z,V(:,k));
%% -------------------------- Finding Charges -----------------------------
tic;
%    cd1
for j=1:cd1
    temp(j,k)=0;
    Q(1,k)= temp(j,k);
    Q(1,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(1,k)=(Q(1,k))/l1;
    fprintf('%d \n',j)
end
fprintf('Q on cd1 Found.\n')

%    cd2
for j=(cd1+1):(cd1+cd2)
    temp(j,k)=0;
    Q(2,k)= temp(j,k);
    Q(2,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(2,k)=(Q(2,k))/l2;
    fprintf('%d \n',j)
end
%    cd3
for j=(cd1+cd2+1):(cd1+cd2+cd3)
    temp(j,k)=0;
    Q(3,k)= temp(j,k);
    Q(3,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(3,k)=(Q(3,k))/h2;
    fprintf('%d \n',j)
end
%    cd4
for j=(cd1+cd2+cd3+1):(cd1+cd2+cd3+cd4)
    temp(j,k)=0;
    Q(4,k)= temp(j,k);
    Q(4,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(4,k)=(Q(4,k))/h2;
    fprintf('%d \n',j)
end
%    cd5
for j=(cd1+cd2+cd3+cd4+1):(cd1+cd2+cd3+cd4+cd5)
    temp(j,k)=0;
    Q(5,k)= temp(j,k);
    Q(5,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(5,k)=(Q(5,k))/h1;
    fprintf('%d \n',j)
end
%    ca1
for j=(cd1+cd2+cd3+cd4+ca1+1):(cd1+cd2+cd3+cd4+cd5+ca1)
    temp(j,k)=0;
    Q(6,k)= temp(j,k);
    Q(6,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(6,k)=(Q(6,k))/l1;
    fprintf('%d \n',j)
end
fprintf('Q on ca1 Found.\n')

%    ca2
for j=(cd1+cd2+cd3+cd4+ca1+ca2+1):Ns
    temp(j,k)=0;
    Q(7,k)= temp(j,k);
    Q(7,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(7,k)=(Q(7,k))/h1;
    fprintf('%d \n',j)
end
%da1
for j=(Ns+1):Ns+da1
    temp(j,k)=0;
    Q(8,k)= temp(j,k);
    Q(8,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(8,k)=(Q(8,k))/l8;
    fprintf('%d \n',j)
end
fprintf('Q on da1 Found.\n')

%da2
for j=(Ns+1)+da1:Ns+da1+da2
    temp(j,k)=0;
    Q(9,k)= temp(j,k);
    Q(9,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(9,k)=(Q(9,k))/l9;
    fprintf('%d \n',j)
end
%da3
for j=(Ns+1)+da1+da2:Ns+da1+da2+da3
    temp(j,k)=0;
    Q(10,k)= temp(j,k);
    Q(10,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(10,k)=(Q(10,k))/a;
    fprintf('%d \n',j)
end
%da4
for j=(Ns+1)+da1+da2+da3:Ns+da1+da2+da3+da4
    temp(j,k)=0;
    Q(11,k)= temp(j,k);
    Q(11,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(11,k)=(Q(11,k))/a;
    fprintf('%d \n',j)
end
%da5
for j=(Ns+1)+da1++da2+da3+da4:Ns+da1+da2+da3+da4+da5
    temp(j,k)=0;
    Q(12,k)= temp(j,k);
    Q(12,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(12,k)=(Q(12,k))/l4;
    fprintf('%d \n',j)
end
%dd1
for j=(Ns+1)+da1+da2+da3+da4+da5:Ns+da1+da2+da3+da4+da5+dd1
    temp(j,k)=0;
    Q(13,k)= temp(j,k);
    Q(13,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(13,k)=(Q(13,k))/b;
    fprintf('%d \n',j)
end
fprintf('Q on dd1 Found.\n')

%dd2
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1:Ns+da1+da2+da3+da4+da5+dd1+dd2
    temp(j,k)=0;
    Q(14,k)= temp(j,k);
    Q(14,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(14,k)=(Q(14,k))/h4;
    fprintf('%d \n',j)
end

%dd3
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2:Ns+da1+da2+da3+da4+da5+dd1+dd2+dd3
    temp(j,k)=0;
    Q(15,k)= temp(j,k);
    Q(15,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(15,k)=(Q(15,k))/h4;
    fprintf('%d \n',j)
end
%dd4
for j=(Ns+1)+da1+da2+da3+da4+da5+dd1+dd2+dd3:N
    temp(j,k)=0;
    Q(16,k)= temp(j,k);
    Q(16,k)= temp(j,k) + (cdensity(j,k))*r(j,5);
    Density_q(16,k)=(Q(16,k))/l5;
    fprintf('%d \n',j)
end
toc
end
fprintf('All Done! \n')