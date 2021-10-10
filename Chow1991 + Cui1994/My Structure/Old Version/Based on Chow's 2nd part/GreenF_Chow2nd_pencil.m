% The function provides the inverse Hankel Transform for the Green's
% Function (including the GPOF method).
% no Plasma
% including the eps0
function G = GreenF_Chow2nd_pencil(xn, zn, xm, zm, dx, dz, epsp, epsr)
syms k
ax = dx/2; az = dz/2;
eps0=8.85*1e-12; epseff = epsr/epsp;
lex= 1*1e-3; lem= 1*1e-3; gap= 0; h=1*1e-3; te = 0.06*1e-3;
N_PM=3; MM=N_PM; Mfix=0; zerolimit=1e-6; L=N_PM;% for The Prony's Method

if (zn < 0)
    fprintf('charges are assumed to be at z > 0! \n')
elseif ( zn >= 0 )
%     if(zm >= zn)
        %% d > h, zm >= d
        % epsrs = epsp;
        K = (1-epseff)/(1+epseff);
        F1 = @(k) (K - exp(-2*h*k))./(1- K * exp(-2*h*k)) - K;
        % Prony method
        x = zeros(1,2*N_PM);
        Ts = 2; p = N_PM;
        for ii = 1: 2*N_PM
           x(1,ii)= double(F1(ii*Ts)); 
        end
        ai = zeros(1,p); bi = zeros(1,p); ri = zeros(1,p); airi = zeros(1,p); %self_matrix1=zeros(1,p);
        % method = 'classic';
        [Amp1, alpha1, freq1, theta1] = matrix_pencil (x,p,Ts);
        for kk = 1:p
            ha = (Amp1(kk)*exp(1j*theta1(kk)))/(exp(Ts*(alpha1(kk) + 1j*2*pi* freq1(kk))));
            hb = Ts*(alpha1(kk) + 1j*2*pi* freq1(kk));
            ai(1,kk) = real(ha) + 1j * imag(ha)
            bi(1,kk) = real(hb) + 1j * imag(hb)
            ri(1,kk) = sqrt( (xm - xn)^2 + (0)^2 + ((bi(1,kk)))^2 );
            airi(1,kk) = ai(1,kk) / ri(1,kk);
%             self_matrix1(1,kk) = (ai(1,kk)) * atanh( a/( sqrt(a^2 + ((bi(1,kk)))^2)));
        end
%         if (n1 ~= m1)
            r10 = sqrt( (xm - xn)^2 + (0)^2 + (zm - zn)^2);
            r10_p = sqrt( (xm - xn)^2 + (0)^2 + (zm + zn)^2);
            if (r10 == 0)
                fprintf('something is wrong');
                G = double((1/(1/(4*pi*eps0))) * (K/r10_p + sum(airi)));
            else
                G = double((1/(4*pi*eps0)) * (1/r10 + K/r10_p + sum(airi)));
            end
            if (isnan(G))
                fprintf('err at n1 = %d, m1 = %d \n', n1,m1)
            end
%         elseif(n1 == m1)
%             Zmm0 = 2*(1 + K)*(2*ax*atan(az/sqrt(ax^2 + az^2)) + az*(-log(-ax + sqrt(ax^2 + az^2)) + log(ax + sqrt(ax^2 + az^2))));
%             G = double((1/(4*pi*(epsp*eps0))) * (Zmm0 + sum(Zmmi)));
%             if (isnan(G))
%                 fprintf('err at n1 = m1 = %d \n', n1)
%             end
%         else
%             fprintf('err! \n')
%         end
end
end
% end