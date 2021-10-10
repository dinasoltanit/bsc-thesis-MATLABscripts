clc; clear; close all;
%% IHT of the GF
% f = 2 kHz, V = 12 kVpp
Vt = 6 * 1e3;
lex= 1*1e-3; lem= 1*1e-3; gap= 0; h=0.17*1e-3; te = 0.06*1e-3;   
%eps0=8.85*1e-12; epsr= 2.7*eps0; epsp = 1*eps0; %epsp=1+77.62*1j;
eps0=8.85*1e-12; epsr= 2.7; epsp = 1; %epsp=1+77.62*1j;
N1 = 10; N2 = 6; N_MOM= N1*N2;
dx=0.1*1e-3; dz=0.01*1e-3; %mm
% if I say, epsp = 1, it means that we are working in air; however, if I
% say that epsp = e1 + 1j * e2, it means that we are working in a plasmonic
% region.
%% Prony's Methos for the sake of GF's coeffs.
% -- first eq.
Z1=zeros(N_MOM,N_MOM); b1=zeros(N_MOM,1); m1=0; n1=0;
b1(:,1) = Vt;
% b1(:,1) = Vt*eps0;
for xm1= 0.05*1e-3 : dx: (lex - 0.05*1e-3)
    for zm1= (h+0.005*1e-3): dz: (h+te - 0.005*1e-3)
        m1 = m1 + 1;
        n1 = 0;
        for xn1= 0.05*1e-3 : dx: (lex - 0.05*1e-3)
            for zn1= (h+0.005*1e-3): dz: (h+te - 0.005*1e-3)
                n1 = n1 + 1;
%                   if (zn1 <= zm1)
                      % 1
                        tic
                        % Z matrix
                        G = Gt2G_up_electrode(zm1, zn1, xm1, xn1, h, epsp, epsr, n1, m1, dx, dz);
                        Z1(m1,n1)= G;
                        toc
%                   elseif ( (h <= zm1)&& (zm1 < zn1) )
%                       % 2
%                         G = Gtilda2G(zm1, zn1, xm1, xn1, h, epsp, epsr, n1, m1, dx, dz);
%                         Z1(m1,n1)= dx*dz*(double(G));
%                   else
%                       fprintf("An error has occured at %d and %d. \n" , m1 , n1)
%                   end
            end
        end
    end   
end
figure(1)
a1= real( (Z1)\b1 );
rho_6kV = reshape((a1),N2,N1);
imagesc(rho_6kV);
colorbar