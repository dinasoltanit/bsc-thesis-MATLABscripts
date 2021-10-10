% The function provides the inverse Hankel Transform for the Green's
% Function (including the GPOF method).
% no Plasma
% including the eps0
function G = GreenF_Chow2nd(xn, zn, xm, zm, dx, dz, epsp, epsr)
syms k
ax = dx/2; az = dz/2;
eps0=8.85*1e-12; epseff = epsr/epsp;
N_PM=3; MM=N_PM; Mfix=1; zerolimit=1e-6; L=N_PM;% for The Prony's Method
h = 0.17 *1e-3;
if (zn < 0)
    fprintf('charges are assumed to be at z > 0! \n')
elseif ( zn >= 0 )
        % epsrs = epsp;
        K = (1-epseff)/(1+epseff)
        F01 = K;
        F1 = @(k) (K - exp(-2*k))./(1- K * exp(-2*k)) - F01;
        % Prony method
        x = zeros(1,2*N_PM);
        y = zeros(1,2*N_PM);
        Ts = 0.5; % defines the equal space between the samples
        for ii = 1: 2*N_PM
            x(1,ii)= (ii) * Ts;
            y(1,ii)= double(F1((ii) * Ts)); 
        end
        [Poles,Coeffs,M,mu] = gpof(y,x,L,MM,Mfix,zerolimit);

            ai = Coeffs
            bi = Poles
            ri = sqrt( (xm - xn)^2 + (0)^2 + (zm + zn - h*(bi)).^2 );
            airi = ai ./ ri;

            r10 = sqrt( (xm - xn)^2 + (0)^2 + (zm - zn)^2);
            r10_p = sqrt( (xm - xn)^2 + (0)^2 + (zm + zn)^2);
            if (r10 == 0)
                fprintf('something is wrong');
                G = double((1/(4*pi*(eps0))) * (K/r10_p + sum(airi)));
            else
                G = double((1/(4*pi*(eps0))) * (1/r10 + K/r10_p + sum(airi)));
            end
            if (isnan(G))
                fprintf('err at n1 = %d, m1 = %d \n', n1,m1)
            end

end
end
% end