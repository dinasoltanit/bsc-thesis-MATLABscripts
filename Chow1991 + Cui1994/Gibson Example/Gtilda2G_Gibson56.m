function G = Gtilda2G_Gibson56(xm, ym, xn, yn, dx, dy, n, m, epsr)
syms x y xp yp
    eps0 = 1;%eps0 = 8.85 * 10^(-12);
    Gs = @(x,y,xp,yp) 1./((( (x - xp).^2 + (y - yp).^2 ).^(0.5)));
%     Gs_int = @(xp,yp) Gs(xm,ym,xp,yp);
%     G_int2 = integral2( Gs_int , (xn-dx/2), (xn+dx/2), (yn-dy/2), (yn+dy/2) );
    if ( n==m )
        a = dx/2;
        G = ( (8*a)/(1) ) * log(1+sqrt(2));
    elseif( n~=m )
        G = double(Gs(xm,ym,xn,yn)) *dx*dy ;
    else
        fprintf('An error has occured. \n')
    end
        
end