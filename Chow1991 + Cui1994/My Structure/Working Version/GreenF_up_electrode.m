% The function provides the inverse Hankel Transform for the Green's
% Function (including the GPOF method).
% no Plasma
% including the eps0
function G = GreenF_up_electrode(xn, zn, xm, zm, dx, dz, epsp, epsr)
syms k
ax = dx/2; az = dz/2;
eps0=8.85*1e-12; h = 0.17*1e-3;
N_PM=3; MM=N_PM; Mfix=1; zerolimit=1e-6; L=N_PM;% for The Prony's Method

if (zn < h)
    fprintf('this function is for the upper electrode. zn has to be greater than h! \n');
elseif ( zn >= h )
    if(zm >= zn)
        %% d > h, zm >= d
        % epsrs = epsp;
        A = @(k) -(epsp + epsr - epsp*exp(2*zn.*k/h) - epsr*exp(2*zn.*k/h) - epsp*exp(2*h.*k/h) + epsr*exp(2*h.*k/h) + epsp*exp(2*k/h.*(zn - h)) - epsr*exp(2*k/h.*(zn - h)))./(4*pi*k/h.*(epsp*exp(2*zn.*k/h) + epsr*exp(2*zn.*k/h) - epsp*exp(2*k/h.*(zn - h)) + epsr*exp(2.*k/h*(zn - h))));
        phi = @(k) (A(k) .* exp(-k/h.*(zm - zn))) *(1/(epsp*eps0))*(-1);
        F0 = double(limit((phi(k) .* 2*(epsp*eps0).*k/h) , k , inf));
        F = @(k) phi(k) .* 2*(epsp*eps0).*k/h - F0;
        % Prony method
        x = zeros(1,2*N_PM);
        y = zeros(1,2*N_PM);
        Ts = 0.5; % defines the equal space between the samples
        for ii = 1: 2*N_PM
            x(1,ii)= (ii) * Ts;
            y(1,ii)= double(F((ii) * Ts)); 
        end
        [Poles,Coeffs,M,mu] = gpof(y,x,L,MM,Mfix,zerolimit);
%         p = M;
%         a1 = zeros(1,p); b1 = zeros(1,p); r1 = zeros(1,p); compims1 = zeros(1,p); self_matrix1=zeros(1,p);
%          disp(mu);
%         for kk = 1:p
            ai = Coeffs;
            bi = Poles;
            ri = sqrt( (xm - xn)^2 + (0)^2 + ((h* bi)).^2 );
            airi = ai ./ ri;
            Zmmi = 2*dz*(ai) .* atanh( ax./( sqrt(ax^2 + ((h* bi)).^2)));
%         end
%         if (n1 ~= m1)
            r10 = sqrt( (xm - xn)^2 + (0)^2 );
            if (r10 == 0)
                G = double((1/(4*pi*(epsp*eps0))) * (sum(airi))) *dx*dz;
            else
                G = double((1/(4*pi*(epsp*eps0))) * (F0/r10 + sum(airi))) *dx*dz;
            end
            if (isnan(G))
                fprintf('err at n1 = %d, m1 = %d \n', n1,m1)
            end
%         elseif(n1 == m1)
%             Zmm0 = dz* ax * F0 * ( log(-ax) + log(ax) )/( sqrt(ax^2) );
%             G = double((1/(4*pi*(epsp*eps0))) * ( Zmm0 + sum(Zmmi) ));
%             if (isnan(G))
%                 fprintf('err at n1 = m1 = %d \n', n1)
%             end
%         else
%             fprintf('err! \n')
%         end

    elseif (zm <= zn && zm >= h)
        %% d > h, zm <= d && zm >= h
        B =@(k) ((exp(-k/h.*(zn-h))).*(epsp*(1-exp(2*h.*k/h))+epsr*(1+exp(2*h.*k/h))))./(4*pi*k/h.*(epsp*(-1+exp(-2*h.*k/h)) + epsr*(1+exp(2*h.*k/h))));
        C =@(k) (exp(-k/h.*(zn-h)))./(4*pi.*k/h);

        phi = @(k) (B(k) .* exp(-k/h.*(zm-h)) + C(k) .* exp(k/h .*(zm-h)))*(1/(epsp*eps0))*(-1);
        F0 = double(limit((phi(k) .* 2*(epsp*eps0).*k/h) , k , inf));
        F = @(k) phi(k) .* 2*(epsp*eps0).*k/h - F0;
        % Prony method
        x = zeros(1,2*N_PM);
        y = zeros(1,2*N_PM);
        Ts = 0.5; % defines the equal space between the samples
        for ii = 1: 2*N_PM
            x(1,ii)= (ii) * Ts;
            y(1,ii)= double(F((ii) * Ts)); 
        end
        [Poles,Coeffs,M,mu] = gpof(y,x,L,MM,Mfix,zerolimit);
%         p = M;
%         a1 = zeros(1,p); b1 = zeros(1,p); r1 = zeros(1,p); compims1 = zeros(1,p); self_matrix1=zeros(1,p);
%          disp(mu);
%         for kk = 1:p
            ai = Coeffs;
            bi = Poles;
            ri = sqrt( (xm - xn)^2 + (0)^2 + ((h* bi)).^2 );
            airi = ai ./ ri;
            Zmmi = 2*dz*(ai) .* atanh( ax./( sqrt(ax^2 + ((h* bi)).^2)));
%         end
%         if (n1 ~= m1)
            r10 = sqrt( (xm - xn)^2 + (0)^2 );
            if (r10 == 0)
                G = double((1/(4*pi*(epsp*eps0))) * (sum(airi))) *dx*dz;
            else
                G = double((1/(4*pi*(epsp*eps0))) * (F0/r10 + sum(airi))) *dx*dz;
            end
            if (isnan(G))
                fprintf('err at n1 = %d, m1 = %d \n', n1,m1)
            end
%         elseif(n1 == m1)
%             Zmm0 = dz* ax * F0 * ( log(-ax) + log(ax) )/( sqrt(ax^2) );
%             G = double((1/(4*pi*(epsp*eps0))) * ( Zmm0 + sum(Zmmi) ));
%             if (isnan(G))
%                 fprintf('err at n1 = m1 = %d \n', n1)
%             end
%         else
%             fprintf('err! \n')
%         end
    elseif (zm <= h && zm >= 0)
        %% d > h, zm <= h && zm >= 0
        D = @(k) (epsp*(-k/h.*(zn - 2*h)))./(2*pi.*k/h .*(epsp*(-1 + exp(2*h.*k/h)) + epsr*(1+exp(2*h.*k/h))));
        
        phi = @(k) (D(k) .* (exp(-k/h.*(zm)) - exp(k/h.*(zm))))*(1/(epsr*eps0))*(-1);
        F0 = double(limit((phi(k) .* 2*(epsp*eps0).*k/h) , k , inf));
        F = @(k) phi(k) .* 2*(epsp*eps0).*k/h - F0;
        % Prony method
        x = zeros(1,2*N_PM);
        y = zeros(1,2*N_PM);
        Ts = 0.5; % defines the equal space between the samples
        for ii = 1: 2*N_PM
            x(1,ii)= (ii) * Ts;
            y(1,ii)= double(F((ii) * Ts)); 
        end
        [Poles,Coeffs,M,mu] = gpof(y,x,L,MM,Mfix,zerolimit);
%         p = M;
%         a1 = zeros(1,p); b1 = zeros(1,p); r1 = zeros(1,p); compims1 = zeros(1,p); self_matrix1=zeros(1,p);
%          disp(mu);
%         for kk = 1:p
            ai = Coeffs;
            bi = Poles;
            ri = sqrt( (xm - xn)^2 + (0)^2 + ((h* bi)).^2 );
            airi = ai ./ ri;
            Zmmi = 2*dz*(ai) .* atanh( ax./( sqrt(ax^2 + ((h* bi)).^2)));
%         end
%         if (n1 ~= m1)
            r10 = sqrt( (xm - xn)^2 + (0)^2 );
            if (r10 == 0)
                G = double((1/(4*pi*(epsp*eps0))) * (sum(airi))) *dx*dz;
            else
                G = double((1/(4*pi*(epsp*eps0))) * (F0/r10 + sum(airi))) *dx*dz;
            end
            if (isnan(G))
                fprintf('err at n1 = %d, m1 = %d \n', n1,m1)
            end
%         elseif(n1 == m1)
%             Zmm0 = dz* ax * F0 * ( log(-ax) + log(ax) )/( sqrt(ax^2) );
%             G = double((1/(4*pi*(epsp*eps0))) * ( Zmm0 + sum(Zmmi) ));
%             if (isnan(G))
%                 fprintf('err at n1 = m1 = %d \n', n1)
%             end
%         else
%             fprintf('err! \n')
%         end
    end
else
    fprintf('errer! \n')
end
end