function [S_approx]=Serii_fb(Fun,a,N)
   syms x;
   a(a == 0) = eps;
   derivative(x)=Fun(x);
   S_approx=subs(derivative,x,a);
   for ii=1:N
      derivative = diff(derivative);
      S_approx = S_approx + subs(derivative,x,a)/factorial(ii);
   end
end