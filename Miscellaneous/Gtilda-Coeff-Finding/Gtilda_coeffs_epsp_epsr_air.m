clc; clear; close all;
% clear A B C D epsp epsr k d h
syms A B C D epsp epsr k d h eps0 W

eq1= B*exp(-k*(d-h))+ C*exp(k*(d-h)) == A;
eq2= (epsp/epsr)*D*(exp(-k*h) - exp(k*h)) - C == B;
eq3= B*exp(-2*k*(d-h)) - (A + (W*eps0 * epsp/k))* exp(-k*(d-h)) == C;
eq4= 2*C * ((epsp/epsr - 1)*exp(-k*h) - (epsp/epsr + 1)*exp(k*h))^(-1) == D;

[solA, solB, solC, solD] = solve([eq1, eq2, eq3, eq4], A, B, C, D);

fA= simplify(solA)
fB= simplify(solB)
fC= simplify(solC)
fD= simplify(solD)
