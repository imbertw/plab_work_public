%dnusolver.m : This code is the discrete resonance dispersion solver that uses sisolver3d2 to find the delta_nu_fsr (GHz), Q, group
%index, FSR (GHZ), mode volume (um^3) of a circular ring.
%function [dFSR_GHz neff t Q groupIndex FSR1_GHz FSR2_GHz Veff_um3] = dnusolver(width, rOut, clambda, xs, nEffGuess)
%TODO: Output beta_fwm
function [dFSR_GHz neff t Q groupIndex FSR1_GHz FSR2_GHz Veff_um3] = dnusolver(width, rOut, clambda, xs, nEffGuess)
tic
c = 299792458;           % m/s
um = 1e-6;               % micron conversion

%Set up vector of wavelengths to be checked

%% Wavelength range

%Estimate how low in wavelength you need to go
blambda = floor(1/(1/clambda+1/(3*pi*rOut))*(1/0.001))/(1/0.001);  %start wavelength rounded to 1nm
%Estimate how high in wavelength you need to go
ulambda = ceil(1/(1/clambda-1/(3*pi*rOut))*(1/0.001))/(1/0.001);   %end wavelength rounded to dlambda
lambdas = blambda:(ulambda-blambda)/4:ulambda;

% Calculate omega in Hz
omega = 2*pi./(lambdas*um)*c;

%% Mode number (azimuthal propagation constant)
gamma = zeros(1,length(lambdas));       %Initialize gamma

%% Set options for sisolver3d2
OPTS.NMODES_CALC = 4;                 %Number of modes you want to calculate
OPTS.eigmode = 'b';                   %beta mode
%PML parameters
OPTS.PMLwidth = [0 0.5 0 0]; % [left right bottom top]
OPTS.PMLsigma = [0.2 0.2];
%Effective index guess
k0 =2*pi/lambdas(1);
OPTS.mu_guess = nEffGuess*k0*rOut;     %guess gamma

OPTS.epsavg = 'simple'; % simple dielectric averaging at step changes
%%
for ii = 1:length(lambdas)
    %Set up index layers
    [nlyrs, dlyrsx, dlyrsy, left_to_rout] = eval([xs '(width,lambdas(ii))']);
    
    rCenter = rOut-width/2.;
    
    %Calculate wavevector
    k0 = 2*pi/lambdas(ii); 
    
    %Initalize modesolver
    dxy           = [0.01, 0.01];                   %Discretization
    
    OPTS.radius   = (rOut-left_to_rout); % specify bent mode simulation and set left side of simulation domain
    %Run modesolver
    [N,F] = sisolver3d4(nlyrs, dlyrsx, dlyrsy, dxy, k0, OPTS);
    
    % Generate SZ for the mode volume, but only for the central frequency
    if ii==3
        N1 = N;
        F1 = F;
        nlyrs1 = nlyrs;
        dlyrsx1 = dlyrsx;
        dlyrsy1 = dlyrsy;
    end
    
    %Find the azimuthal propagation constant gamma
    gamma(1,ii)     = F.beta(1,1);
    OPTS.mu_guess   = real(F.beta(1,1));
    %Convert from angular propagation constant to linear propagation
    %constant
    beta(ii) = gamma(1,ii)/rCenter;       % beta is the longitudinal propagation constant

end

%% Interpolate the angular frequencies for finer sampling
%resample the lines
dg = 0.01; 
xgamma = (min(real(gamma(1,:)))+1):dg:max(real(gamma(1,:))-1);
pomega = spline(real(gamma), real(omega),xgamma);                          %Pump frequencies
somega = spline(real(gamma+1), real(omega),xgamma);                        %Signal frequencies (where there is a pump, signal, and idler)
iomega = spline(real(gamma-1), real(omega),xgamma);                        %Idler frequencies (where there is a pump, signal, and idler)

%% Calculate FSR difference vs gamma
dnu      = (2*pomega-somega-iomega)/(2*pi*1e9);     %FSR mismatches in GHz
plambda  = 2*pi./pomega*c*1e9;                      %Pump wavelength in nm
dFSR_GHz = spline(plambda,dnu,clambda*1e3);         %FSR mismatch at 1550

%% Calculate the FSR
fsr1     = (pomega-somega)/(2*pi*1e9);       %FSRs in GHz
fsr2     = (iomega-pomega)/(2*pi*1e9); 
FSR1_GHz = spline(plambda,fsr1,clambda*1e3); %FSRs at 1550
FSR2_GHz = spline(plambda,fsr2,clambda*1e3);

%% Calculate the effective index
neff = real(beta(3))/(2*pi/clambda);   %Effective index at center wavelength
%% Calculate the mode volume of the ring
%[Veff_um3]        = mode_volume_ring(SZ, width, rOut); % The mode volume in microns^3
n_nl = index_Si_CMG(clambda); % index of silicon for this frequency, based on which version of you're using
G = [F1, F1, F1, F1]; %just use the same mode, as signal pump and idler modes are all pretty similar
nmodes = [1 1 1 1];   % use the fundamental mode
%[Veff_um3]   = nonlinear_mode_volume(N1, G, nmodes, region_nl, n_nl, clambda); %the nonlinear interaction volume in microns^3
[Veff_um3]   = nonlinear_mode_volume(N1, G, nmodes, n_nl, clambda); %the nonlinear interaction volume in microns^3
%% Calculate the group index of the ring
%resample omega
dOmega       = (max(omega)-min(omega))/100;
xOmega_plus  = omega(3) + dOmega;
xOmega_minus = omega(3) - dOmega;
%Find the gradient
beta_plus  = spline(omega,real(beta),xOmega_plus);
beta_minus = spline(omega,real(beta),xOmega_minus);

betas = [beta_minus real(beta(3)) beta_plus];

dBetadOmega      =  diff(betas*1e6)/dOmega;    %s/m

groupVelocity    = 1./dBetadOmega;             %m/s
avgGroupVelocity = mean(groupVelocity);
groupIndex       = c/avgGroupVelocity;            % group index n_g

%% Calculate the loss Q of the ring
alpha = imag(beta(3));
Q = k0*groupIndex/(2*alpha);          %Q is the quality factor

%% end time
t = toc;
end