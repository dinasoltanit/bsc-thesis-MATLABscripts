function [Amp, alpha, freq, theta] = matrix_pencil (xt,p,Ts)
x = xt';
N = length(x);
Y = hankel( x(1:end-p) , x(end-p:end) );

Y1 = Y(:,1:end-1);
Y2 = Y(:,2:end);

l = eig(pinv(Y1) * Y2); % eigenvalues

alpha = log( (abs(l)) /Ts);
freq = atan2(imag(l) , real(l))/(2*pi*Ts);

Z = zeros(N,p);
for ii = 1:length(l)
   Z(:,ii) = transpose(l(ii) .^(0:N-1)); 
end
rZ = real(Z);
iZ = imag(Z);
% here, Inf values are substituted by realmax values
rZ(isinf(rZ)) = realmax * sign(rZ(isinf(rZ)));
iZ(isinf(iZ)) = realmax * sign(iZ(isinf(iZ)));

Z = rZ + 1i * iZ;
h = Z\x;
Amp = abs(h);
theta = atan2(imag(h) , real(h));
end