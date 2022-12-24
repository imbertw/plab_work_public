%---------.---------.---------.---------.---------.---------.---------.--------|
% Code:      2D PML Vector-Field Modesolver
% Component: PML Matrix Generator using Complex Coordinate Stretching
% Date:      June 7, 2002
% Author:    Milos Popovic
%
% [Pr, Pz, Rmtx] = m2dpmlmatx(w, er(2*M-1,2*N-1), Rvec(2*M-1), Zvec(2*N-1), ...
%                  [sigmaRMAX sigmaZMAX], [dPMLRlo dPMLRhi dPMLZlo dPMLZhi]);
%
% (All input parameters are required)
%
% Field grid is M x N (Er is (M-1)xN, Ez is Mx(N-1)).
% Dielectric er is (2M-1)x(2N-1) as are Pr, Pz and the radius matrix, Rmtx.

% Jun  7, 2002 - First written code.
% Jun 27, 2002 - Correction: radius matrix Rmtx augmented to contain the complex
%                r-coordinate within the PML region! (without this, it was only
%                a "perfect PML" at large radii!). The complex component is
%                in dRmtx, and is returned within Rmtx; such that the total
%                complex radial variable Rmtx = Rmtxold + dRmtx/(omega*epsilon).
%---------.---------.---------.---------.---------.---------.---------.--------|

function [Pr, Pz, Rmtx] = m2dpmlmatx(w, er, Rvec, Zvec, sigmaMAX, dPML)
% Physical constants
c = 299792458; u0 = 4e-7 * pi; e0 = 1/c^2/u0;     % Permittivity, permeability.

% PML thickness at inner(1)/outer(2) radial and lower(3)/upper(4) z- boundaries
dPMLRlo = dPML(1); dPMLRhi = dPML(2); dPMLZlo = dPML(3); dPMLZhi = dPML(4);

% Total and useful (without PML) computational domain edges.
RminT = min(Rvec); RmaxT = max(Rvec); ZminT = min(Zvec); ZmaxT = max(Zvec);
Rmin = RminT + dPMLRlo; Rmax = RmaxT - dPMLRhi;   % Radial PML setup
Zmin = ZminT + dPMLZlo; Zmax = ZmaxT - dPMLZhi;   % Vertical PML setup

% PML max conductivities for parabolic conductivity profile
sigmaRMAX = sigmaMAX(1);      % Max conductivity in radial PML layers
sigmaZMAX = sigmaMAX(2);      % Max conductivity in vertical PML layers

% Complex coordinate stretching (PML) factor
fR = i * sigma(Rvec, [Rmin Rmax], sigmaRMAX) * ...
         ones(size(Zvec(:).')) ./ (w*e0*er);      % Radial sigma at (m,n)
fZ = i * ones(size(Rvec(:))) * ...                % Vertical sigma at (m,n)
         sigma(Zvec, [Zmin Zmax], sigmaZMAX).' ./ (w*e0*er);

% PML Matrices
Pr   = 1 ./ (1 + fR);
Pz   = 1 ./ (1 + fZ);
Rmtx = Rvec(:) * ones(size(Zvec(:).'));
% [MP] Jun 27/02: Adding dRmtx complex coord; Rtotal = Rmtx + dRmtx/(omega*eps)
dRmtx= rcomplex(Rvec, [Rmin Rmax], sigmaRMAX) * ones(size(Zvec(:).'));
Rmtx = Rmtx + (dRmtx ./ (w*e0*er));     % Total complex rho coordinate (Jun27).


%---------.---------.---------.---------.---------.---------.---------.--------|
% FUNCTION s(:) = sigma(...): Parabolic conductivity function for PML layer    |
%---------.---------.---------.---------.---------.---------.---------.--------|
function s = sigma(Rvec, Rlimits, sigmaMAX) %, dPML)
Rmin = Rlimits(1); Rmax = Rlimits(2);             % PML-"free domain" boundaries
dPMLlo = Rmin - min(Rvec); dPMLhi = max(Rvec) - Rmax;   % Find PML thicknesses

% Set up conductivity profile
s = zeros(size( Rvec(:) ));
iL = find(Rvec < Rmin);
if (dPMLlo > 0) & ~isempty(iL)
    sL = sigmaMAX * ((Rmin-Rvec)/dPMLlo).^2;
    s(iL) = sL(iL);
end
iR = find(Rvec > Rmax);
if (dPMLhi > 0) & ~isempty(iR)
    sR = sigmaMAX * ((Rvec-Rmax)/dPMLhi).^2;
    s(iR) = sR(iR);
end

%-[MP] Jun 27, 2002-.---------.---------.---------.---------.---------.--------|
% FUNCTION r(:) = rcomplex(...): Complex r-coord *component* within PML layer  |
%---------.---------.---------.---------.---------.---------.---------.--------|
function r = rcomplex(Rvec, Rlimits, sigmaMAX) %, dPML)
Rmin = Rlimits(1); Rmax = Rlimits(2);             % PML-"free domain" boundaries
dPMLlo = Rmin - min(Rvec); dPMLhi = max(Rvec) - Rmax;   % Find PML thicknesses

% Set up regular coordinate
r = zeros(size(Rvec));
iL = find(Rvec < Rmin);             % Complex-stretch the coordinate within PML
if (dPMLlo > 0) & ~isempty(iL)
    rL = i * sigmaMAX/3 * (Rvec-Rmin).^3/(dPMLlo.^2); % Missing r + and 1/(w*eps)
    r(iL) = rL(iL);
end
iR = find(Rvec > Rmax);
if (dPMLhi > 0) & ~isempty(iR)
    rR = i * sigmaMAX/3 * (Rvec-Rmax).^3/(dPMLhi.^2); % Missing r + and 1/(w*eps)
    r(iR) = rR(iR);
end
r = r(:);                           % Output a column vector.
