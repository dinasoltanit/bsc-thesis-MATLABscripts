function [Poles,Coeffs,M,mu]=gpof(y,x,L,MM,Mfix,zerolimit)
% Author: Amir Ahmad Shishegar
% Date: Tir 5, 1380.
% Version 1.0
% Status: Works fine.
%===============================================================================
% help:
%
% GPOF method to find coeffs and poles of exponential approximation expansion:
% y(x)=sum k=1 ^M [coeffs(k)*exp(poles(k)*x)]
%
% Input parameters:
% y   :  A vector containing N samples calculated at equally spaced points.
% x   :  A vector with N equally spaced elements that y is calculated in.
%                       (x(i)-x(i-1)=constant) x and y can be complex.
% L   :  An integer that can be changed to make better approximation. It should
%        satisfies  M<=L<=(N-M) to have the condition of the main theorem.
%                       Some tests show that L=N/2 is the best approximation.
% MM  :  An integer that can be set to show the number of exponents required in
%        the approximation. It is effective only if Mfix is set to 1.
% Mfix:  An integer. If Mfix=1, the M sets to MM else M calculated based on the
%                       zerolimit and variable found in the program.(automatically)
% zerolimit : A positive real number that is used to select M. It is effective
%                       if and only if Mfix~=1. 1e-6 seems good for that.
%
% Output parameters:
% poles  : The poles of the system in a Mx1 vector.
% coeffs : The coefficients of the exponents in a Mx1 vector.
% M      : The order of the system (number of poles) calculated if Mfix~=1.
% mu     : The exponential approximation Nx1 vector calculated in x locations.
%          If the approximation is done well, it must be close to y.
%
%  This program is based on Hua's and Sarkar's method.
%===========================Routin begins here.=================================
x=x(:);
deltat=x(2)-x(1);
y=y(:);
[m,n]=size(y);
if n~=1
    error('Input must be a vector');
end
N=m;
%=================== Generating Y1 and Y2 from input data vector y. =====
Y1=zeros(N-L,L);
ii=1:N-L;
for jj=1:L,
    Y1(:,jj)=y(ii+jj-1);
end;

Y2=Y1(:,2:L);
Y2=[Y2 y(L+1:N)];
%=================== Finding the poles of the system. ==========================
% 1- Singular value decomposition of Y1.
[U,S,V]=svd(Y1);
% 2- Decision making about M. If Mfix=1 then M=MM
% else, M calculated based on algorithm below.
if min(size(S))~=1
    Sdiag=diag(S);
else
    Sdiag=S(1);
end
zerovalue=zerolimit*Sdiag(1);
if Mfix==1
    M=MM;
else
    temp=(1:min(N-L,L))';
    temp=(Sdiag>zerovalue).*temp;
    [~,M]=max(temp); % M is known now.
end
% 2'- Checking if M,N,L have the required condition in the theorem or not.
if (L<M)||(L>(N-M))
    error('N, L and M should satisfy relation: M<=L<=N-M')
end
% 3- Dimension changing of U, S and V versus M.
if M==0
    error('M=0');
end
Sinv=diag(1./Sdiag(1:M));    % new S : a MxM matrix is inversed.
Un=U(:,1:M);                 % new U : a (N-L)xM matrix.
Vn=V(:,1:M);                 % new V : a LxM matrix.
% 4- Ze calculation.
Ze=Sinv*(Un'*Y2*Vn);
% 5- Finding the eigenvalues of Ze.
try
    zi=eig(Ze);                  % zi is a vector containing M eigenvalues.
catch err
    if strcmp(err.identifier,'MATLAB:eig:matrixWithNaNInf')
        zi = ones(M,1);
    else
        rethrow(err)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important:
% zi can be calculated using generalized eigenvalue calculation of matlab:
% Y1p=pinv(Y1); % Pseudoinverse of Y1.
% A=Y1p*Y2;
% B=Y1p*Y1;
% zzi=eig(A,B); % This is a vector of length L containing all generalized
%               % eigenvalues of A and B.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6- Calculating the poles.
Poles=log(zi)./deltat;
%=================== Finding the coefficients. =================================
% 1- z1 making. It is a (N-L)xM matrix.
z1=zeros(N-L,1);
z1(1:M)=zi;
z1=fliplr(vander(z1)).';
z1=z1(:,1:M);
% 2- z2 making. It is a MxL matrix.
z2=zeros(L,1);
z2(1:M)=zi;
z2=fliplr(vander(z2));
z2=z2(1:M,:);
% 3- Coefficients calculation.
Coeffs=pinv(z1)*Y1*pinv(z2);
Coeffs=diag(Coeffs);
% 4- Scaling of coeffs to find the coefficients of the main approximation.
Coeffs=Coeffs./exp(Poles*x(1));
%=================== Test of the results.=======================================
mu=exp(x*Poles.')*Coeffs;

%============== sort results based on real and imaginary parts of poles ===
temp = Poles;
temp(abs(real(Poles))<1e-3*abs(Poles)) = 1i*imag(Poles(abs(real(Poles))<1e-3*abs(Poles)));
temp(abs(imag(Poles))<1e-3*abs(Poles)) = real(Poles(abs(imag(Poles))<1e-3*abs(Poles)));
[~,IX] = sort(temp);
Poles = Poles(IX);
Coeffs = Coeffs(IX);
return