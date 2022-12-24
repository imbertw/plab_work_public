% pair_rate.m : This is a function that allows you to
% determine the nonlinear performance of your microring. It calculates:
% The photon flux in units of photons per second based on Cale's equation
% in p.84 of his thesis, given a pump power of 1mW.
%function [I_coinc] = pair_rate(clambda, dnu, Q, Veff)
% inputs: clambda in microns, dnu in GHz, Veff in um^3.
% output in pairs/second
% assumes Kerr index from Cale's optica paper, and assumes we're using Si
% in 45RFSOI

function [I_coinc] = pair_rate(clambda, dnu, Q, Veff)
c =299792458;
um = 1e-6;
c_um = c/um; %units of microns per second
n2 = 3.96e-6; %units of microns^3
P = 1e-3;
frequency = c / (clambda*um);   % in Hz
omega = 2*pi*frequency;         % in radians/s
dnu_Hz = dnu*1e9;               % in Hz
domega = 2*pi*dnu_Hz;           % in radians/s
I_coinc = 1.4622*omega*n2^2*c_um^2*Q*P^2./(index_Si_CMG(clambda)^2*Veff.^2.*(domega.^2+49/9*omega.^2./Q.^2));

