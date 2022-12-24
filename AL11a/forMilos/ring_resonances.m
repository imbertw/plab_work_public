%ring_resonances.m - given a ring geometry (width, radius) and a center
%frequency, this code will find the actual resonant wavelengths around that
%center frequency.
%function [omegas gammas SZ] = ring_resonances(width, rOut, clambda, xs, nEffGuess)
function [lambda_res omegas gammas SZ] = ring_resonances(width, rOut, clambda, xs, nEffGuess)
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
OPTS.PMLwidth = [0 0.5 0.5 0.5]; % [left right bottom top]
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
    dxy           = [0.02, 0.015];                   %Discretization
    
    OPTS.radius   = (rOut-left_to_rout); % specify bent mode simulation and set left side of simulation domain
    %Run modesolver
    [N,F] = sisolver3d2(nlyrs, dlyrsx, dlyrsy, dxy, k0, OPTS);
    
    % Generate SZ for the mode volume, but only for the central frequency
%     if ii==3
%         SZ.N = N;
%         SZ.F = F;
%     end
    SZ(ii).N = N;
    SZ(ii).F = F;
    %modeview(SZ)
    
    %Find the azimuthal propagation constant gamma
    gamma(1,ii)     = F.beta(1,1);
    OPTS.mu_guess   = real(F.beta(1,1));
    %Convert from angular propagation constant to linear propagation
    %constant
    beta(ii) = gamma(1,ii)/rCenter;       % beta is the longitudinal propagation constant

end

%% Interpolate the angular frequencies for finer sampling
%resample the lines
%dg = 0.01; 
%xgamma = (min(real(gamma(1,:)))):dg:max(real(gamma(1,:)));
%pomega = spline(real(gamma), real(omega),xgamma);                          %Pump frequencies

%% Interpolate to obtain the angular resonances
%Find resonant mode numbers (whole numbers)
round_min_gamma = ceil(min(real(gamma(1,:))));
round_max_gamma = floor(max(real(gamma(1,:))));
gammas = round_min_gamma:1:round_max_gamma;
omegas = spline(real(gamma),real(omega),gammas);
lambda_res = 2*pi./(omegas*um)*c;


y = ones(1,length(lambda_res));
figure;
stem(lambda_res,y,'LineWidth',2);

%% end time
t = toc
end
