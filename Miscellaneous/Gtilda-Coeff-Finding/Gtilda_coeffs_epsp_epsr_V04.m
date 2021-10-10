clc; clear; close all;
% clear A B C D epsp epsr k d h
syms A B C D epsp epsr k d h eps0 W

eq1= B*exp(-k*(h-d))+ C*exp(k*(h-d)) == A;
eq2= D*(exp(-k*d) - exp(k*d)) - C == B;
eq3= ((-1)*epsp/epsr)*A*exp(-k*(h-d)) + B*exp(-2*k*(h-d)) == C;
eq4= ((-1/(4*pi*k))- C )*exp(-k*d) == D;

[solA, solB, solC, solD] = solve([eq1, eq2, eq3, eq4], A, B, C, D);

fA= simplify(solA)
fB= simplify(solB)
fC= simplify(solC)
fD= simplify(solD)
