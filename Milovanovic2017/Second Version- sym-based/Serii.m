function [S_approx]=Serii(Fun,a,N,dd,zmm)
   syms x;
   d = dd;
   zm = zmm;
   a(a == 0) = eps;
   derivative(x)=Fun(x,d,zm);
   S_approx=subs(derivative,x,a);
   for ii=1:N
      derivative = diff(derivative);
      S_approx = S_approx + subs(derivative,x,a)/factorial(ii);
   end
end