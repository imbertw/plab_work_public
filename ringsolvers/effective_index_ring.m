%neff_check.m : this code uses sisolver3d2 to obtain the effective index
%of a structure given the cross-section xs, the radius, and the width.
%% Start of function
%function [neff SZ] = effective_index(width, rOut, clambda, xs, nEffGuess)
function [neff SZ] = effective_index_ring(width, rOut, clambda, xs, nEffGuess)
tic;
c = 299792458;           % m/s
lambda0 = clambda;          % microns
[nlyrs, dlyrsx, dlyrsy, left_to_rout] = eval([xs '(width,lambda0)']);
k0 = 2*pi/lambda0; 
rCenter = rOut-width/2.;
%Run modesolver
dxy = [0.02, 0.005];                   %Discretization
OPTS.NMODES_CALC = 4;                 %Number of modes you want to calculate
OPTS.eigmode = 'b';                   %beta mode
OPTS.PMLwidth = [0 0.5 0 0]; % [left right bottom top]
OPTS.PMLsigma = [0.2 0.2];
OPTS.epsavg = 'simple'; % simple dielectric averaging at step changes
OPTS.mu_guess = nEffGuess*k0*rOut;     %guess gamma
if rOut ~= inf;
    OPTS.radius = (rOut-left_to_rout); % specify bent mode simulation and set left side of simulation domain
end

[N,F] = sisolver3d2(nlyrs, dlyrsx, dlyrsy, dxy, k0, OPTS);
%make imaginary values of N.y, N.x and N.n have color 0
Ny = N.y;
Nx = N.x;
Nn = N.n;
Ny(imag(Ny) ~= 0) = 0;
Nx(imag(Nx) ~= 0) = 0;
Nn(imag(Nn) ~= 0) = 0;
imagesc(Nx, Ny, (Nn).'); set(gca, 'YDir', 'normal'); colorbar; % Take a look at the structure
%Convert from angular propagation constant to linear propagation
%constant
% gamma = real(F.beta(1,1)); %change the row index to pick out the correct mode
% beta = gamma/rCenter;

gamma= real(F.beta(1:OPTS.NMODES_CALC,1));
beta = gamma/rCenter;
% figure; imagesc(N.x, N.y, (N.n).'); set(gca, 'YDir', 'normal');
% colorbar; %creates a picture of the layer stack

SZ.F = F; SZ.N = N;
neff = beta/(2*pi/lambda0);   %Effective index at center wavelength
beta = sort(beta,1,'descend'); %re-order 
neff = sort(neff,1,'descend');
neff = neff(1);
toc;




 

 