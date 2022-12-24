function [beta_fwm] = betafwm(clambda, V);
n = index_Si_CMG(clambda);
n2 = 3.96*1e-14*1e-4; %in units of m^2/W
c = 299792458; % in units of ms^-1
V = V*1e-18; %mode volume in m^3
beta_fwm = n2*c./(n^2*V);