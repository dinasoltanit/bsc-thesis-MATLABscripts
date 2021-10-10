clc;clear;
%-----------------------Parameters--------------------------------------------------------------------------------------------------
% z_prime=10;
% h1=-10;
% h2=-20;
% eps1=1.2;
% eps2=1.5;
% syms z k
syms k h1 h2 z z_prime eps1 eps2 r J
%-----------------------Defining the solution: v_tilda------------------------------------------------------------------------------
fA =(eps1*exp(2*k*(h1 + h2)) + eps2*exp(2*k*(h1 + h2)) - eps1*exp(2*k*(h1 + z_prime)) + eps1*exp(2*k*(h2 + z_prime)) - eps2*exp(2*k*(h1 + z_prime)) - eps2*exp(2*k*(h2 + z_prime)) - eps1*exp(2*h1*k) + eps1*exp(2*h2*k) - eps2*exp(2*h1*k) - eps1*exp(4*h1*k) - eps2*exp(2*h2*k) + eps2*exp(4*h1*k) - eps1^2*exp(2*k*(h1 + h2)) - eps1*exp(2*k*(2*h1 + z_prime)) + eps2*exp(2*k*(2*h1 + z_prime)) + eps1^2*exp(2*k*(h1 + z_prime)) - eps1^2*exp(2*k*(h2 + z_prime)) + eps1*exp(2*k*(h1 + h2 + z_prime)) + eps2*exp(2*k*(h1 + h2 + z_prime)) - eps1^2*exp(2*h1*k) + eps1^2*exp(2*h2*k) + eps1^2*exp(4*h1*k) - eps1^2*exp(2*k*(2*h1 + z_prime)) + eps1^2*exp(2*k*(h1 + h2 + z_prime)) + eps1*eps2*exp(2*k*(h1 + h2 + z_prime)) - eps1*eps2*exp(2*k*(h1 + h2)) + eps1*eps2*exp(2*k*(h1 + z_prime)) + eps1*eps2*exp(2*k*(h2 + z_prime)) - eps1*eps2*exp(2*h1*k) - eps1*eps2*exp(2*h2*k) - eps1*eps2*exp(4*h1*k) + eps1*eps2*exp(2*k*(2*h1 + z_prime)))/(4*k*pi*(eps1*exp(2*k*(h2 + z_prime)) - eps1*exp(2*k*(h1 + z_prime)) - eps2*exp(2*k*(h1 + z_prime)) - eps2*exp(2*k*(h2 + z_prime)) - eps1*exp(2*k*(2*h1 + z_prime)) + eps2*exp(2*k*(2*h1 + z_prime)) + eps1^2*exp(2*k*(h1 + z_prime)) - eps1^2*exp(2*k*(h2 + z_prime)) + eps1*exp(2*k*(h1 + h2 + z_prime)) + eps2*exp(2*k*(h1 + h2 + z_prime)) - eps1^2*exp(2*k*(2*h1 + z_prime)) + eps1^2*exp(2*k*(h1 + h2 + z_prime)) + eps1*eps2*exp(2*k*(h1 + h2 + z_prime)) + eps1*eps2*exp(2*k*(h1 + z_prime)) + eps1*eps2*exp(2*k*(h2 + z_prime)) + eps1*eps2*exp(2*k*(2*h1 + z_prime))));
fB1 = -(eps1*exp(2*h1*k) - eps2*exp(2*k*(h1 + h2)) - eps1*exp(2*k*(h1 + h2)) - eps1*exp(2*h2*k) + eps2*exp(2*h1*k) + eps1*exp(4*h1*k) + eps2*exp(2*h2*k) - eps2*exp(4*h1*k) + eps1^2*exp(2*k*(h1 + h2)) + eps1^2*exp(2*h1*k) - eps1^2*exp(2*h2*k) - eps1^2*exp(4*h1*k) + eps1*eps2*exp(2*k*(h1 + h2)) + eps1*eps2*exp(2*h1*k) + eps1*eps2*exp(2*h2*k) + eps1*eps2*exp(4*h1*k))/(4*k*pi*(eps1*exp(k*(2*h1 + 2*h2 + z_prime)) + eps2*exp(k*(2*h1 + 2*h2 + z_prime)) - eps1*exp(k*(2*h1 + z_prime)) + eps1*exp(k*(2*h2 + z_prime)) - eps2*exp(k*(2*h1 + z_prime)) - eps1*exp(k*(4*h1 + z_prime)) - eps2*exp(k*(2*h2 + z_prime)) + eps2*exp(k*(4*h1 + z_prime)) + eps1^2*exp(k*(2*h1 + 2*h2 + z_prime)) + eps1^2*exp(k*(2*h1 + z_prime)) - eps1^2*exp(k*(2*h2 + z_prime)) - eps1^2*exp(k*(4*h1 + z_prime)) + eps1*eps2*exp(k*(2*h1 + 2*h2 + z_prime)) + eps1*eps2*exp(k*(2*h1 + z_prime)) + eps1*eps2*exp(k*(2*h2 + z_prime)) + eps1*eps2*exp(k*(4*h1 + z_prime))));
fB2 = exp(-k*z_prime)/(4*k*pi); 
fC1 =-(eps1*exp(3*h1*k) + eps2*exp(3*h1*k) - eps1*exp(k*(h1 + 2*h2)) + eps2*exp(k*(h1 + 2*h2)))/(2*k*pi*(eps1*exp(k*(2*h1 + 2*h2 + z_prime)) + eps2*exp(k*(2*h1 + 2*h2 + z_prime)) - eps1*exp(k*(2*h1 + z_prime)) + eps1*exp(k*(2*h2 + z_prime)) - eps2*exp(k*(2*h1 + z_prime)) - eps1*exp(k*(4*h1 + z_prime)) - eps2*exp(k*(2*h2 + z_prime)) + eps2*exp(k*(4*h1 + z_prime)) + eps1^2*exp(k*(2*h1 + 2*h2 + z_prime)) + eps1^2*exp(k*(2*h1 + z_prime)) - eps1^2*exp(k*(2*h2 + z_prime)) - eps1^2*exp(k*(4*h1 + z_prime)) + eps1*eps2*exp(k*(2*h1 + 2*h2 + z_prime)) + eps1*eps2*exp(k*(2*h1 + z_prime)) + eps1*eps2*exp(k*(2*h2 + z_prime)) + eps1*eps2*exp(k*(4*h1 + z_prime))));
fC2 =(eps2*exp(3*h1*k) - eps1*exp(3*h1*k) + eps1*exp(k*(h1 + 2*h2)) + eps2*exp(k*(h1 + 2*h2)))/(2*k*pi*(eps1*exp(k*(2*h1 + 2*h2 + z_prime)) + eps2*exp(k*(2*h1 + 2*h2 + z_prime)) - eps1*exp(k*(2*h1 + z_prime)) + eps1*exp(k*(2*h2 + z_prime)) - eps2*exp(k*(2*h1 + z_prime)) - eps1*exp(k*(4*h1 + z_prime)) - eps2*exp(k*(2*h2 + z_prime)) + eps2*exp(k*(4*h1 + z_prime)) + eps1^2*exp(k*(2*h1 + 2*h2 + z_prime)) + eps1^2*exp(k*(2*h1 + z_prime)) - eps1^2*exp(k*(2*h2 + z_prime)) - eps1^2*exp(k*(4*h1 + z_prime)) + eps1*eps2*exp(k*(2*h1 + 2*h2 + z_prime)) + eps1*eps2*exp(k*(2*h1 + z_prime)) + eps1*eps2*exp(k*(2*h2 + z_prime)) + eps1*eps2*exp(k*(4*h1 + z_prime)))); 
fD = -(eps1*exp(k*(2*h1 + h2)))/(k*pi*(eps1*exp(k*(2*h1 + 2*h2 + z_prime)) + eps2*exp(k*(2*h1 + 2*h2 + z_prime)) - eps1*exp(k*(2*h1 + z_prime)) + eps1*exp(k*(2*h2 + z_prime)) - eps2*exp(k*(2*h1 + z_prime)) - eps1*exp(k*(4*h1 + z_prime)) - eps2*exp(k*(2*h2 + z_prime)) + eps2*exp(k*(4*h1 + z_prime)) + eps1^2*exp(k*(2*h1 + 2*h2 + z_prime)) + eps1^2*exp(k*(2*h1 + z_prime)) - eps1^2*exp(k*(2*h2 + z_prime)) - eps1^2*exp(k*(4*h1 + z_prime)) + eps1*eps2*exp(k*(2*h1 + 2*h2 + z_prime)) + eps1*eps2*exp(k*(2*h1 + z_prime)) + eps1*eps2*exp(k*(2*h2 + z_prime)) + eps1*eps2*exp(k*(4*h1 + z_prime))));
A=((fA));
B1=((fB1));
B2=((fB2));
C1=((fC1));
C2=((fC2));
D=((fD));
%-----------------------------------------------------------------------------------------------------------------------------------
H = piecewise(z>z_prime,A*exp(-k*(z-z_prime)), (0<z)&(z<z_prime), B1*exp(-k*z)+B2*exp(k*z), -h1<z<0, C1*exp(-k*(z+h1))+C2*exp(k*(z+h1)), -h2<z<-h1, D*(exp(-k*(z+h2))+exp(k*(z+h2))))

%---------------------------Inverse Hankel Transform--------------------------------------------------------------------------------
[G,J]=iht(H,k,r,J)
%---------------------------Integral Equation---------------------------------------------------------------------------------------
% syms x b
% solve(int(f(x,b),x,t0,t1)-L,b)
%------------------------------------------------------------------------------------------------------------------------------------
%Inverse Hankel transform of order 0.
%Input:
% H      Spectrum K(k)
% k      Spatial frequencies [rad/m]   {pi/numel(H)*(0:numel(H)-1)}
% r      Radial positions [m]          {0:numel(H)-1}
% J      Integration kernel °          {default}
%
%Output:
% h      Signal h(r)
% J      Integration kernel
%------------------
function [G,J]=iht(H,k,r,J)
if sum(size(H) > 1) > 1
   error('Spectrum must be a vector.');
end
if nargin < 2 || isempty(k)
   k=pi/numel(H)*(0:numel(H)-1).';
else
   [k,w]=sort(k(:).');
   H=H(w);
end
if nargin < 3 || isempty(r)
   r=0:numel(H)-1;
end
if nargin < 4 || isempty(J)
   k=[(k(2:end) + k(1:end-1))/2 k(end)];
   J=1/2/pi./r(:)*k.*besselj(1,r(:)*k);
   J(r == 0,:)=1/4/pi*k.*k;
   J=J - [zeros(numel(r),1) J(:,1:end-1)];
elseif exist('w','var')
   J=J(:,w);
end
G=reshape(J*H(:),size(r));
end

