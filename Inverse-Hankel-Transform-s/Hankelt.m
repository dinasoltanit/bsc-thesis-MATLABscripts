function H=Hankelt(a,b,r,dt,n)
n2=2*n;
%J0_x=besselj(0,x);
%J1_x=besselj(1,x);
d=(a+b)/(2*b) + 1e-10; % to avoid 0^0 in u(x,k)
% u_xk= (x-d)^(k-1);
% u_xk_prime= (k-1)*(x-d)^(k-2);
%f_x=1/(x^2+1);
g_x=0;
point= a+ ((1:dt:n) - 1)*(b-a)/(n-1);
rhs=zeros(1,n2);
for i=1:dt:n
    rhs(1,i)=point(i)*f_x(point(i));
end
for i=1:dt:n
   rhs(1,n+i)=g_x; 
end
mat=zeros(n2,n2);
for j=1:dt:n
    for k=1:dt:n
        mat(j,k)=u_xk_prime(point(j),k,d);
    end  
end
for j=1:dt:n
    for k=n+1:dt:n2
        mat(j,k)=r*u_xk(point(j),k-n,d);
    end  
end
for j=1:dt:n
    for k=1:dt:n
        mat(j+n,k)=-r*u_xk(point(j),k,d);
    end  
end
for j=1:dt:n
    for k=n+1:dt:n2
        mat(j+n,k)=u_xk_prime(point(j),k-n,d)-u_xk(point(j),k-n,d)/point(j);
    end  
end
c=rhs/mat;
for i=1:dt:n
   sum1=sum(c(i)*u_xk(b,i,d));
   sum2=sum(c(i)*u_xk(a,i,d));
   sum3=sum(c(n+i)*u_xk(b,i,d));
   sum4=sum(c(n+i)*u_xk(a,i,d));
end
H=(sum1*besselj(0,r*b) - sum2*besselj(0,r*a) + sum3*besselj(1,r*b)- sum4*besselj(1,r*a))/n;

end