clc; clear; close all;
epsr = 2.7; % for Kapton
f0 = 2*1e3; w0 = 2*pi*f0; e= 1.6*1e-19; eps0 = 8.85*1e-12;
P = 1; nu_m = 5.3e9*P*760; m_e = 9.1094e-31; n = 1e19;
wp = sqrt((n*e^2)/(e0 * m_e));
epsp = 1-((wp/f0)^2)/(1+(nu_m/f0)*1j);
eps = [epsr epsp 1];
for ii = 1:3
    
    for jj = 1:3
        
    
    
    end
end




















