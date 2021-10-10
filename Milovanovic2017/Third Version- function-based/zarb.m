function [aj]=zarb(a,jj)
if jj==0
    aj=1;
elseif jj>0
   for n=1:jj
       aj=aj*(a-n+1);
    end
end    
end