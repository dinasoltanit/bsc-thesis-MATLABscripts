function [Amp, alpha, freq, theta] = polynomial_method (xt,p,Ts, method)
% method: 'classic', 'ls' or 'tls' (case sensitive)
%% define the solving method
    CLASSIC = 0;
    LS = 1;
    TLS = 2;
    x = xt';
    N = length(x);
    
    %fprintf('N = %d', N);
    %disp(x);
    if strcmpi (method, 'classic')
        if N ~= 2*p
            disp ('ERROR: Length of x must be 2*p samples in classical method.');
            Amp = []; alpha = []; freq = []; theta = [];
            return;
        else
            solve_method = CLASSIC;
        end
    elseif strcmpi (method, 'ls')
        solve_method = LS;
    elseif strcmpi (method, 'tls')
        solve_method = TLS;
    else
        disp ('ERROR: error in parsing the argument "method".');
        Amp = []; alpha = []; freq = []; theta = [];
        return;
    end
    %% step 1
    T = toeplitz( x(p:N-1) , x(p:-1:1) );
    switch solve_method
        case {CLASSIC , LS}
            x_vect = (x(p+1:N));
            a = - T \ x_vect;
        case TLS
            a = tls(T, -x(p+1:N));
    end
    %% check for indeterminate forms
    indeterminate_form = (sum(isnan(a)) || sum(isinf(a)));
    if (indeterminate_form)
        Amp = []; alpha = []; freq = []; theta = [];
        return;
    end
    
    %% step 2
    c = transpose([1;a]);
    r = roots(c);
    alpha = log(abs(r))/Ts;
    freq = atan2(imag(r) , real(r))/(2*pi*Ts);
    % In case alpha equals to +/-Inf the signal will not be recovered for
    % n=0, since Inf*0 = NaN. Making alpha = +/-realmax that
    % interditerminance will be solved.
    alpha(isinf(alpha)) = realmax * sign(alpha(isinf(alpha)));
    
    %% step 3
    switch solve_method
        case CLASSIC
            len_vandermonde = p; % exact case (N=2p) find h with p samples
        case LS
            len_vandermonde = N; % overdetermined case (N>2p) find h with N samples
        case TLS
            len_vandermonde = N; % overdetermined case (N>2p) find h with N samples
    end
    
    Z = zeros(len_vandermonde, p);
    
    for i = 1:length(r)
        Z(:,i) = transpose( r(i) .^(0:(len_vandermonde-1)) );
    end
    
    rZ = real(Z);
    iZ = imag(Z);
    % here, Inf values are substituted by realmax values
    rZ(isinf(rZ)) = realmax * sign(rZ(isinf(rZ)));
    iZ(isinf(iZ)) = realmax * sign(iZ(isinf(iZ)));
    
    Z = rZ + 1i * iZ;
    
    switch solve_method
        case {CLASSIC, LS}
            h = Z\(x(1:len_vandermonde));
        case TLS
            % if there exists NaN values, SVD won't work
            indeterminate_form = sum( sum( isnan(Z) || isinf(Z) ) );
            if (indeterminate_form)
                Amp = []; alpha = []; freq = []; theta = [];
                return;
            else
                h = tls(Z,x(1:len_vondermonde));
            end
    end
    Amp = abs(h);
    theta = atan2(imag(h) , real(h));
end