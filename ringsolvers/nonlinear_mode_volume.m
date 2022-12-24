%nonlinear_mode_volume_si.m: function that takes SZ from modesolver and
%returns the nonlinear interaction volume in cubic microns
%Also returns the beta_fwm in units of Joules^-1 as defined in Savvy's fwm_memo2.pdf
%This particular code assumes the material is silicon

%inputs:
%SZ: output of mode solver
%width: width of ring in microns
%rOut: outer radius of ring in microns
%n_nl: linear refractive index of the nonlinear material

%outputs:
%Veff_um3: nonlinear interaction volume in cubic microns

function [Veff_um3] = nonlinear_mode_volume(N, F, nmodes, n_nl, wavelength)

% calculate the effective 3rd order nonlinear effect coefficients: \alpha,
% \beta, which both have unit of 1/energy
% for example, da3/dt = ... -iw3* beta3* a1^3. 
% inputs: 
%         F, the four field profiles of eigenmode of linear operators
%         N, refractive index distribution matrix
%         region_nl = [xleft, xright, yleft, yright], the region of
%                      nonlinear material
%         chi3 = [chi3_1111, chi3_1122, chi3_1212, chi3_1221]; the four independent elements for Si at
%                 wavelength far away from resonance
%         nmmodes = [n1, n2, n3, n4], which modes to use in F; 
% 06/27/2011: for [100] wafer; and 
n_2 = 3.96e-18;             %Kerr index of silicon in m^2/W
c = 299792458;              %speed of light in m/s
chi3ratio = [0.4233];     % chi3_1122/chi3_1111, % for silicon near 1.55um; the value 0.4233 comees from Lin's paper; 
epsilon_0 = 8.85e-12; 

% old version with 3 B's; 
% [A, B1122, B1212, B1221] = AB_nl3w(F, region_nl, [1 1 1 1], 'R', Lcavity); 
% norm_energy = sqrt(prod([energyfcn(N, F(1), nmodes(1)), energyfcn(N, F(2), nmodes(2)), energyfcn(N, F(3), nmodes(3)), energyfcn(N, F(4), nmodes(4))]));
% coeff = 1/16* epsilon_0* (chi3a * A + chi3b*B1122 +chi3c*B1212 +chi3d*B1221)/ norm_energy;   % in unit of per joule; 

[A, B] = AB_nl3w4(F, N, n_nl, nmodes, 'R');   % Savvy 062711: updated AB coefficient taking into account of the fact that field direction is angle dependent for ring mode

norm_energy = sqrt(prod([energyfcn(N, F(1), nmodes(1)), energyfcn(N, F(2), nmodes(2)), energyfcn(N, F(3), nmodes(3)), energyfcn(N, F(4), nmodes(4))]));

Veff_um3 = 4*norm_energy/(epsilon_0^2*n_nl^4*(A+B*chi3ratio));
%beta_fwm = 1e18*n_2*c/(index_Si_CMG(wavelength)^2*Veff_um3); %1e18 is for the conversion of mode volume from cubic microns to cubic meters
