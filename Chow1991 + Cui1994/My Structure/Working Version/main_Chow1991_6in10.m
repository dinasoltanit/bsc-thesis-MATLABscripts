clc; clear; close all;
%% IHT of the GF
% f = 2 kHz, V = 12 kVpp
Vt = 6 * 1e3;
lex= 1*1e-3; lem= 1*1e-3; gap= 0; h=0.17*1e-3; te = 0.06*1e-3;   
% eps0=8.85*1e-12; epsr= 2.7*eps0; epsp = 1*eps0; %epsp=1+77.62*1j;
eps0=8.85*1e-12; epsr= 2.7; epsp = 1; %epsp=1+77.62*1j;
N1 = 10; N2 = 6; N_MOM= N1*N2;
dx= 2*0.1*1e-3; dz=3*0.01*1e-3; %mm


xx = 0:dx:2*1e-3;
zz = h:dz:h+1.5*1e-3;
Voltage = zeros(numel(xx),numel(zz));
[X,Z] = meshgrid(xx,zz);
mode = 'noPlasma';
 m = 0; n = 0;
for xxx = 0:dx:2*1e-3
    m = m + 1;
    n = 0;
    for zzz = h:dz:h+1.5*1e-3
        n = n + 1;
        tic
        V = voltage_calc(xxx,zzz,mode)
        Voltage(m,n) = V;
        toc
    end
end
Ex = zeros(numel(xx),numel(zz)); Ez = zeros(numel(xx),numel(zz));
for m = 1:numel(xx)
    for n = 1:numel(zz)
        if m == 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if n == 1
                Ex(m,n)= ( -Voltage(m+2,n) + 4*Voltage(m+1,n) - 3*Voltage(m,n) )/(2*dx);
                Ez(m,n)= ( -Voltage(m,n+2) + 4*Voltage(m,n+1) - 3*Voltage(m,n) )/(2*dz);
            elseif n== numel(zz)
                Ex(m,n)= ( -Voltage(m+2,n) + 4*Voltage(m+1,n) - 3*Voltage(m,n) )/(2*dx);
                Ez(m,n)= ( Voltage(m,n-2) - 4*Voltage(m,n-1) + 3*Voltage(m,n) )/(2*dz);
            else
                Ex(m,n)= ( -Voltage(m+2,n) + 4*Voltage(m+1,n) + -3*Voltage(m,n) )/(2*dx);
                Ez(m,n)= ( Voltage(m,n+1) - Voltage(m,n-1) )/(2*dz);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif m== numel(xx)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if n == 1
                Ex(m,n)= ( Voltage(m-2,n) - 4*Voltage(m-1,n) + 3*Voltage(m,n) )/(2*dx);
                Ez(m,n)= ( -Voltage(m,n+2) + 4*Voltage(m,n+1) - 3*Voltage(m,n) )/(2*dz);
            elseif n== numel(zz)
                Ex(m,n)= ( Voltage(m-2,n) - 4*Voltage(m-1,n) + 3*Voltage(m,n) )/(2*dx);
                Ez(m,n)= ( Voltage(m,n-2) - 4*Voltage(m,n-1) + 3*Voltage(m,n) )/(2*dz);
            else
                Ex(m,n)= ( Voltage(m-2,n) - 4*Voltage(m-1,n) + 3*Voltage(m,n) )/(2*dx);
                Ez(m,n)= ( Voltage(m,n+1) - Voltage(m,n-1) )/(2*dz);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if n == 1
                Ex(m,n)= ( Voltage(m+1,n) - Voltage(m-1,n) )/(2*dx);
                Ez(m,n)= ( -Voltage(m,n+2) + 4*Voltage(m,n+1) - 3*Voltage(m,n) )/(2*dz);
            elseif n== numel(zz)
                Ex(m,n)= ( Voltage(m+1,n) - Voltage(m-1,n) )/(2*dx);
                Ez(m,n)= ( Voltage(m,n-2) - 4*Voltage(m,n-1) + 3*Voltage(m,n) )/(2*dz);
            else
                Ex(m,n)= ( Voltage(m+1,n) - Voltage(m-1,n) )/(2*dx);
                Ez(m,n)= ( Voltage(m,n+1) - Voltage(m,n-1) )/(2*dz);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

figure(1)
pcolor(X,Z,Ex)
title('Ex distribution (V/m)')
xlabel('along z-direction (m)')
ylabel('along x-direction (m)')
figure(2)
pcolor(X,Z,Ez)
title('Ez distribution (V/m)')
xlabel('along z-direction (m)')
ylabel('along x-direction (m)')
figure(3)
pcolor(X,Z,sqrt(Ex.^2 + Ez.^2))
title('E magnitude distribution (V/m)')
xlabel('along z-direction (m)')
ylabel('along x-direction (m)')
figure(4)
pcolor(X,Z,Voltage)
title('Voltage distribution (V)')
xlabel('along z-direction (m)')
ylabel('along x-direction (m)')