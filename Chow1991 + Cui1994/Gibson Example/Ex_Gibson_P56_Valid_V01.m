clc; clear; close all;
%% IHT of the GF
% f = 2 kHz, V = 12 kVpp
Vt = 1;
L = 1; epsr = 1; eps0=8.85*10^(-12); % eps0 = 1;
N1=10; N2=10;
N_MOM= N1*N2;
dx=L/N1; dy=L/N2; %mm
%% Prony's Methos for the sake of GF's coeffs.
% -- first eq.
Z1=zeros(N_MOM,N_MOM); b1=zeros(N_MOM,1); m=0; n=0;
b1(:,1) = 4*pi*epsr*eps0;
for xm = L/(2*N1):dx:(1-L/(2*N1))
    for ym = L/(2*N2):dy:(1-L/(2*N2))
        m = m + 1;
        n = 0;
        for xn = L/(2*N1):dx:(1-L/(2*N1))
            for yn = L/(2*N2):dy:(1-L/(2*N2))
                n = n + 1;
            tic
            % Z matrix
            Z1(m,n) = Gtilda2G_Gibson56(xm, ym, xn, yn, dx, dy, n, m, epsr);
            toc
            end
        end
    end
end
a1= ( Z1\b1 ) * 10^12;
a = reshape(a1,N1,N2);
imagesc(a);
colorbar