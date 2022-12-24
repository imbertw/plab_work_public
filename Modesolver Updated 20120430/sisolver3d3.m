% Step-index 2D-cross-section waveguide modesolver bootstrap function
% Milos Popovic, Mar 21, 2005
%
% Syntax:   [N,F,V] = sisolver3d(nlyrs, dlyrsx, dlyrsy, dxy, k0, OPTS)
%
% Inputs:   nlyrs           [NxM] matrix of section refractive indices along y (1:N) and x (1:M)
%
%           dlyrsx          [1xM] vector of M section widths along x
%           dlyrsy          [1xN] vector of N section heights along y
%           dxy             [1x2] [dx dy] discretizations along x and y
%           k0              Free space k-vector (2*pi/free space wavelength)
%           (all 4 of the above spatial parameters must use the same units)
%
%           OPTS            [OPTIONAL] Any options to be passed to solver engine (e.g. m2wcyl)
%                           Format:  OPTS = struct(name, value, name, value, ...);
%                           Example: OPTS = struct('epsavg','pn','radius',5, ...);
%
%               OPTS.epsavg     Chooses 'pn' (default, pn = parallel/normal field)
%                               or 'simple' averaging of the dielectric constant
%                               at step interfaces.
%               OPTS.radius     Set radius of left interface (of entire
%                               computational domain) if bend mode is to be computed.
%                               NOTE: leftmost radial coordinate must be >0.
%               OPTS.hsymm      Set to 'E' or 'M' for electric or magnetic
%                               wall horizontal symmetry (i.e. about a vertical axis)
%               OPTS.vsymm      Set to 'E' or 'M' for electric or magnetic
%                               wall vertical symmetry (i.e. about a horizontal axis)
%               OPTS.NMODES_CALC [Default = 1] How many modes to calculate
%               OPTS.NMODES_KEEP [1 to NMODES_CALC; default = NMODES_CALC] The first how many
%                               modes to return as result
%               OPTS.enginever  [Default = 'm2wcylR2'] allows calling fcns to be
%                               told which version of solver engine to call (old = 'm2wcyl',
%                               newer = 'm2wcylR2')
%
%               OPTS.BC        = Default = [0 0 0 0]; [left right bottom top] boundary conditions,
%                                 PEC = 0; PMCh = 1; PBC = 2; PMC = 3; PECh = 4.
%                                 PEC  = perfect electric conductor, defined at pixel edge
%                                 PMC  = perfect magnetic conductor, defined at pixel edge
%                                 PECh = perfect electric conductor, defined at pixel center
%                                 PMCh = perfect magnetic conductor, defined at pixel center
%                                 PBC  = periodic boundary condition
%                               Alternative specification of BC's: {'PEC','PEC','PMCh','PMC'}
%
%               OPTS.<other>:   any other solver engine-accepted OPTS options can
%                               also be passed in.  See help of solver engine for its
%                               OPTS parameters.
%
% Outputs:
%           N               Refractive index distribution simulated.  For
%                           MxN pixel grid, with dx x dy pixels,
%                             N.x (1:2M+1) x-coordinate vector, every dx/2
%                             N.y (1:2N+1) y-coordinate vector, every dy/2
%                             N.n (2M+1 x 2N+1) refractive index on dx/2,dy/2 grid
%                             (different E field components sample different parts, 1/4 is unused)
%           F               Computed mode field components and coordinates
%           V               Raw eigenvector matrix
%
% Debugging:
%           global VMODE; VMODE = 2;  % provides full debugging info while the solver is running
%
% Examples:
%
% Straight waveguide:
% [N,F] = sisolver3d([1 1 1; 1 3.5 1; 1 1 1], [1 0.6 1], [1 0.3 1], [.02 .02], 2*pi/1.55);  % For 1 mode, or for 3 of them plus a guess index for speed:
% [N,F] = sisolver3d([1 1 1; 1 3.5 1; 1 1 1], [1 0.6 1], [1 0.3 1], [.02 .02], 2*pi/1.55, struct('NMODES_CALC',3,'mu_guess',2.5 * 2*pi/1.55));
% Correct eff. index result:  neff = F.beta/(2*pi/1.55) = [2.7723 2.2624 1.7310];
% (due to coarse grid above this is slightly off the exact result)
%
% Bent waveguide:
% [N,F] = sisolver3d([1 1 1; 1 3.48 1; 1 1 1], [0.7 0.48 1], [1 0.22 1], [.02 .02], 2*pi/1.55, struct('radius',1-0.48/2-0.7,'mu_guess',3*2*pi/1.55*1,'PMLwidth',[0 0.5 0 0],'PMLsigma',[0.2 0.2]));
%
% Note: to view output modes using "modeview", type:
%     modeview( struct('N',N,'F',F) );
% or
%     SZ.F = F; SZ.N = N; modeview(SZ)

%         \|/
%         @ @     Kilroy wuz here!
% ---oOOO-(_)-OOOo-----------------
%
%
% TO DO:        - Make rounding of domain size to discretization more sensible/intuitive?
%               - Remove muguess, NMODES_CALC and NMODES_KEEP as variables,
%                 leave them in the OPTS structure and don't copy them at
%                 start of function.
%
% Code updates:
% -------------
% Mar 21, 2005  - First written.
% Nov 22, 2005  - Bug fix: occasional appearance of NaN in generated index eliminated.
%                 Solution was to replace isinf(npsqy) with ~isfinite(npsqy), since the
%                 purposefully used Inf at edges, later removed from the matrix, also
%                 generate some NaN's when Inf/Inf is encountered.
%                 Problem was: index has NaN if use command
%                    [N]=sisolver3d([1.446*[1 1 1]; 1.446 3.476 1.446; 1.446*[1 1 1]], ...
%                                   [1 .45 1], [1 .22 1], [.02 .02], 2*pi/1.55);
% Jul 19, 2006  - Added allowed arbitrary PML parameters to be passed in
%                 (no longer hard coded).
% Jul 26, 2006  - Can pass in boundary conditions in OPTS.BC.
% Aug 20, 2006  - Added warning when computing bends that radius is for
%                 left edge of computational domain.
% Oct 25, 2006  - [unfinished] Adding support for Star-P
% Oct 30, 2006  - IN PROGRESS: started adding vertical/horizontal domain symmetries
%                 (electric/magnetic) - ported from sisolver3d_symm.m (now
%                 the sisolver3d_symm.m function can be retired).
% Nov  3, 2006  - Made fieldmode, coordmode, eigmode, etc. have default
%                 rather than hard-coded values.  This was done to remove
%                 remaining hard-coded hacks in this code.  Still need to
%                 unhardcode the PMLwidth and PMLsigma
% Nov  4, 2006  - Added (optional) support for new version of solver engine,
%                 m2wcylR2.m, through the OPTS.enginever option.
% Dec 17, 2010  - Added example in text only.
% Feb  4, 2011  - Added computation of bend loss in post-processing
% Apr 28, 2011  - Cleaning up code based on feedback from Tymon Barwicz
%               - Renamed file to sisolver3d2.m as temporary new branch -
%                 many changes will be made, so until we verify no bugs
%                 were introduced by the changes, this will be a separate
%                 branch.
%               - Modified help text at top (clarified index line, deleted
%                 OPTS.PARALLEL option as it wasn't being used, clarified
%                 help for OPTS)
%               - Reorganized some variables, solver engine call, and added
%                 comments
%               - Cleaned up code to make it at Matlab R2011b standards


function [N,F,V] = sisolver3d3(nlyrs, dlyrsx, dlyrsy, dxy, k0, OPTS)
global VMODE; if isempty(VMODE), VMODE = 0; end     % Check diagnostic flag VMODE
if (nargin < 6), OPTS = struct; end;                % Default, just make a blank OPTS structure
if ~isfield(OPTS,'enginever'), OPTS.enginever = 'm2wcylR3'; end             % Default modesolver engine to use
if ~isfield(OPTS,'fieldmode'), OPTS.fieldmode = 'V'; end                    % Default fieldmode is 'V'ectorial (rather than semivectorial 'MX' or 'MY', or scalar 'S')
if ~isfield(OPTS,'coordmode'), OPTS.coordmode = 'C'; end                    % Default coordmode is 'C'artesian (rather than 'R'adial)
if ~isfield(OPTS,'eigmode'),   OPTS.eigmode   = 'b'; end                    % Default eigmode   is 'b'eta propagation constant (rather than 'w' = omega, resonance frequency)
if ~isfield(OPTS,'BC'),        OPTS.BC = [0 0 0 0];  end
if ~isfield(OPTS,'PMLwidth'),  OPTS.PMLwidth = [0 0 0 0];  end              % [default value units: um] [MP20060719] If passed in, keep passed in values, else assign defaults:
if ~isfield(OPTS,'PMLsigma'),  OPTS.PMLsigma = [1 1];      end              % [default value units: 1/um]
if(isfield(OPTS,'mu_guess')),     mu_guess    = OPTS.mu_guess;     else  mu_guess    = max(nlyrs(:)) * k0;  end  % Guess eigenvalue?
if(isfield(OPTS,'NMODES_CALC')),  NMODES_CALC = OPTS.NMODES_CALC;  else  NMODES_CALC = 1;  end                   % Number of modes to find?
if(isfield(OPTS,'NMODES_KEEP')),  NMODES_KEEP = OPTS.NMODES_KEEP;  else  NMODES_KEEP = NMODES_CALC;  end         % Number of modes to keep
% OPTS.eigsfcn = 'eigsmp'; OPTS.v = 40;                                     % Choice of Matlab eigensolver: eigsmp = Milos modified; eigs = original
% OPTS.eigsfcn = 'eigs'; OPTS.v = 40;                                       % Default Matlab solver (v = # of Arnoldi vectors during iterations)
if(VMODE>0), fprintf('2006-04-13 temporary mod: Using eigs instead of eigsmp (because of changed hard drive).\n'); end
% OPTS.eigsfcn = 'jdqr'; %OPTS.v = 40;                                      % Jacobi-Davidson solver

% Parse simpler boundary condition specification method [Apr 30, 2012]
if iscell(OPTS.BC)
    tmp = [0 0 0 0];
    for kk = 1:4
        switch OPTS.BC{kk}
            case 'PEC'
                tmp(kk) = 0;
            case 'PMCh'
                tmp(kk) = 1;
            case 'PBC'
                tmp(kk) = 2;
            case 'PMC'
                tmp(kk) = 3;
            case 'PECh'
                tmp(kk) = 4;
            otherwise
                error('Unrecognized boundary condition specified.');
        end
    end
    OPTS.BC = tmp; clear tmp;
end

% Parse the input structure
%[MM,NN] = size(nlyrs);                                                    % Step-index regions
dlyrsx = dlyrsx(:).'; dlyrsy = dlyrsy(:).';
MM = length(dlyrsx); NN = length(dlyrsy);
dx = dxy(1); dy = dxy(2);

dlyrsx(end) = round(sum(dlyrsx)/dx)*dx + 0*1e-3*dx - sum(dlyrsx(1:end-1));
dlyrsy(end) = round(sum(dlyrsy)/dy)*dy + 0*1e-3*dy - sum(dlyrsy(1:end-1)); % [MP-TODO] Here: round xint(end), yint(end) domain *endpoints only* to be even in grid spacing or else round discretization
xint = [0 cumsum(dlyrsx)]; yint = [0 cumsum(dlyrsy)];       % Interface coordinates
if(VMODE > 0)
    fprintf('Domain limits: (x,y) = user passed (0..%f, 0..%f), rounding last layer to pixel grid (0..%f, 0..%f) | Pixel grid size = %d x %d\n', sum(dlyrsx), sum(dlyrsy), xint(end), yint(end), round([xint(end)/dx yint(end)/dy]));
    fprintf('   x-interface positions at: ');   for k = 1:length(xint), fprintf('%f  ', xint(k)); end;  fprintf('\n');
    fprintf('   y-interface positions at: ');   for k = 1:length(yint), fprintf('%f  ', yint(k)); end;  fprintf('\n');
end
x = (0 : dx/2 : xint(end)); y = (0 : dy/2 : yint(end));     % Grid coordinates
xint(1) = -inf; xint(end) = +inf; yint(1) = -inf; yint(end) = +inf;     % [MP] Put edge interfaces to infinity to avoid index averaging at domain edges

if ((nargin > 5) && isfield(OPTS,'epsavg'))
    if strcmp(OPTS.epsavg,'simple')
        flagepsavg = 0;
    else
        flagepsavg = 1;         % Use parallel/normal field determined epsilon averaging for smooth distribution
    end
else
    flagepsavg = 1;
end

if(flagepsavg == 0)   % Standard (area-arithmetic) index averaging for pixels, then those are averaged for pixel edge tensor components for Ex, Ey, Ez
if(VMODE>0), fprintf('Using simple arithmetic dielectric averaging.\n'); end
% Generate index distribution with area-weighted dielectric averaging (justified by Ampere's law) at interfaces
xp = x(2:2:end-1); yp = y(2:2:end-1); npsq = zeros(length(xp), length(yp));
for mm = 1:MM                                                                           % along x-coordinate
    ixin = find((xp+dx/2) > xint(mm) & (xp-dx/2) <= xint(mm+1));                        % Find all pixels with at least some part inside region (mm,nn)
    xpmin = max(xint(mm), xp(ixin)-dx/2);  xpmax = min(xint(mm+1), xp(ixin)+dx/2);      % Min and max coords of pixel *in* region (mm,nn)
    for nn = 1:NN                                                                       % along y-coordinate
        iyin = find((yp+dy/2) > yint(nn) & (yp-dy/2) <= yint(nn+1));
        ypmin = max(yint(nn), yp(iyin)-dy/2);  ypmax = min(yint(nn+1), yp(iyin)+dy/2);  % Min and max coords of pixel *in* region (mm,nn)

        npsq(ixin,iyin) = npsq(ixin,iyin) + (xpmax-xpmin).'*(ypmax-ypmin)/(dx*dy) * nlyrs(nn,mm)^2;     % Fill index pixel-matrix with area-averaging (near interfaces)
    end
end

% Two-grid index matrix and interface averaging
N.x = x; N.y = y; N.n = zeros(length(x), length(y));
N.n(2:2:end,2:2:end) = sqrt(npsq);                                              % Fill denser two-grid index matrix
N.n(3:2:end-2,:) = sqrt( (N.n(2:2:end-3,:).^2 + N.n(4:2:end-1,:).^2)/2 );       % Dielectric averaging between two pixels at interfaces
N.n(:,3:2:end-2) = sqrt( (N.n(:,2:2:end-3).^2 + N.n(:,4:2:end-1).^2)/2 );
N.n(1,:) = N.n(2,:); N.n(end,:) = N.n(end-1,:); N.n(:,1) = N.n(:,2); N.n(:,end) = N.n(:,end-1);
end

if(flagepsavg == 1)   % More sophisticated index averaging - separate for each tensor component, arithmetic for parallel field, geometric for normal
if(VMODE > 0), fprintf('Using "parallel-normal field" determined dielectric averaging.\n'); end
N.x = x; N.y = y; N.n = zeros(length(x), length(y));
% Better averaging (hopefully at most only 2 pixels shaded near any one interface)
% Ez dielectric const averaging (x -> arithmetic, y -> arithmetic)
xp = x(1:2:end); yp = y(1:2:end); npsqz = zeros(length(xp), length(yp));
for mm = 1:MM                                                                           % along x-coordinate
    ixin = find((xp+dx/2) > xint(mm) & (xp-dx/2) <= xint(mm+1));                        % Find all pixels with at least some part inside region (mm,nn)
    xpmin = max(xint(mm), xp(ixin)-dx/2);  xpmax = min(xint(mm+1), xp(ixin)+dx/2);      % Min and max coords of pixel *in* region (mm,nn)
    for nn = 1:NN                                                                       % along y-coordinate
        iyin = find((yp+dy/2) > yint(nn) & (yp-dy/2) <= yint(nn+1));
        ypmin = max(yint(nn), yp(iyin)-dy/2);  ypmax = min(yint(nn+1), yp(iyin)+dy/2);  % Min and max coords of pixel *in* region (mm,nn)

        npsqz(ixin,iyin) = npsqz(ixin,iyin) + (xpmax-xpmin).'*(ypmax-ypmin)/(dx*dy) * nlyrs(nn,mm)^2;     % Fill index pixel-matrix with area-averaging (near interfaces)
    end
end
%warning off MATLAB:divideByZero
% Ey dielectric const averaging (first y -> geometric, then x -> arithmetic)
xp = x(1:2:end); yp = y(2:2:end-1); npsqy = Inf * ones(length(xp), length(yp)); npsqytot = zeros(length(xp), length(yp));
for mm = 1:MM                                                                           % along x-coordinate
    ixin = find((xp+dx/2) > xint(mm) & (xp-dx/2) <= xint(mm+1));                        % Find all pixels with at least some part inside region (mm,nn)
    xpmin = max(xint(mm), xp(ixin)-dx/2);  xpmax = min(xint(mm+1), xp(ixin)+dx/2);      % Min and max coords of pixel *in* region (mm,nn)
    for nn = 1:NN                                                                       % along y-coordinate
        iyin = find((yp+dy/2) > yint(nn) & (yp-dy/2) <= yint(nn+1));
        ypmin = max(yint(nn), yp(iyin)-dy/2);  ypmax = min(yint(nn+1), yp(iyin)+dy/2);  % Min and max coords of pixel *in* region (mm,nn)

        npsqy(ixin,iyin) = 1./(1./npsqy(ixin,iyin) + 1./((xpmax-xpmin).'*(1./(ypmax-ypmin))/(dx/dy) * nlyrs(nn,mm)^2));     % Fill index pixel-matrix with area-averaging (near interfaces)
    end
    npsqy(~isfinite(npsqy)) = 0;
    npsqytot = npsqytot + npsqy;     % Fill index pixel-matrix with area-averaging (near interfaces)
    npsqy = Inf * ones(length(xp), length(yp));
end
% Ex dielectric const averaging (first x -> geometric, then y -> arithmetic)
xp = x(2:2:end-1); yp = y(1:2:end); npsqx = Inf * ones(length(xp), length(yp)); npsqxtot = zeros(length(xp), length(yp));
for nn = 1:NN                                                                       % along y-coordinate
    iyin = find((yp+dy/2) > yint(nn) & (yp-dy/2) <= yint(nn+1));
    ypmin = max(yint(nn), yp(iyin)-dy/2);  ypmax = min(yint(nn+1), yp(iyin)+dy/2);  % Min and max coords of pixel *in* region (mm,nn)
    for mm = 1:MM                                                                           % along x-coordinate
        ixin = find((xp+dx/2) > xint(mm) & (xp-dx/2) <= xint(mm+1));                        % Find all pixels with at least some part inside region (mm,nn)
        xpmin = max(xint(mm), xp(ixin)-dx/2);  xpmax = min(xint(mm+1), xp(ixin)+dx/2);      % Min and max coords of pixel *in* region (mm,nn)

        npsqx(ixin,iyin) = 1./(1./npsqx(ixin,iyin) + 1./(1./(xpmax-xpmin).'*((ypmax-ypmin))/(dy/dx) * nlyrs(nn,mm)^2));     % Fill index pixel-matrix with area-averaging (near interfaces)
    end
    npsqx(~isfinite(npsqx)) = 0;
    npsqxtot = npsqxtot + npsqx;     % Fill index pixel-matrix with area-averaging (near interfaces)
    npsqx = Inf * ones(length(xp), length(yp));
end
%warning on MATLAB:divideByZero
% Combine tensor components into full index distribution:
N.n(1:2:end,1:2:end) = sqrt(npsqz); N.n(1:2:end,2:2:end-1) = sqrt(npsqytot); N.n(2:2:end-1,1:2:end) = sqrt(npsqxtot);
% Index components on below line should really not matter! and should be NaN, but they *seem* to be used, perhaps in PMLs (coord stretch) so I set them here
N.n(2:2:end-1,2:2:end-1) = sqrt((N.n(1:2:end-2,2:2:end-1).^2 + N.n(3:2:end,2:2:end-1).^2 + N.n(2:2:end-1,1:2:end-2).^2 + N.n(2:2:end-1,3:2:end).^2)/4);
end

if(find(~isfinite(N.n))), error('sisolver3d error: Generated index distribution has non-finite elements.  Report problem to milos@mit.edu.'); end


if(nargout > 1)     % If modes are also desired (F output parameter asked), then call modesolver

% Process horizontal or vertical electric or magnetic wall symmetries (OPTS.vsymm = 'E','M')
if(isfield(OPTS,'vsymm'))   % Asking for any vertical symmetries?
    if(VMODE>0), fprintf('Symmetry condition specified is...: '); end
    switch upper(OPTS.vsymm)
        case 'M'            % Vertical magnetic wall at y = 0 (USER MUST ENSURE MIDDLE LAYER HAS ODD # OF PIXELS!)
            if(VMODE>0), fprintf('Magnetic wall (at half pixel steps) vertical symmetry.\n'); end
            iy = find(N.y > max(N.y)/2 -(dy/2 * 1.001));
            N.y = N.y(iy); N.n = N.n(:,iy);
            OPTS.BC(3) = 1;                 % Set magnetic wall at bottom
        case 'E'            % Vertical electric wall at y = 0 (USER MUST ENSURE MIDDLE LAYER HAS EVEN # OF PIXELS!)
            if(VMODE>0), fprintf('Electric wall (at full pixel steps) vertical symmetry.\n'); end
            iy = find(N.y > max(N.y)/2 -(dy/2 * 0.001));
            N.y = N.y(iy); N.n = N.n(:,iy);
            OPTS.BC(3) = 0;                 % Set electric wall at bottom
        otherwise
            if(VMODE>0), fprintf('unknown (''%s''). ERROR.\n', OPTS.vsymm); end
            error('[MP] sisolver3d - unknown symmetry condition specified');
    end
end
if(isfield(OPTS,'hsymm'))   % Asking for any horizontal symmetries?
    if(VMODE>0), fprintf('Symmetry condition specified is...: '); end
    switch upper(OPTS.hsymm)
        case 'M'            % Horizontal magnetic wall at x = 0 (USER MUST ENSURE MIDDLE LAYER HAS ODD # OF PIXELS!)
            if(VMODE>0), fprintf('Magnetic wall (at half pixel steps) horizontal symmetry.\n'); end
            ix = find(N.x > max(N.x)/2 -(dx/2 * 1.001));
            N.x = N.x(ix); N.n = N.n(ix,:);
            OPTS.BC(1) = 1;                 % Set magnetic wall at left
        case 'E'            % Horizontal electric wall at x = 0 (USER MUST ENSURE MIDDLE LAYER HAS EVEN # OF PIXELS!)
            if(VMODE>0), fprintf('Electric wall (at full pixel steps) horizontal symmetry.\n'); end
            ix = find(N.x > max(N.x)/2 -(dx/2 * 0.001));
            N.x = N.x(ix); N.n = N.n(ix,:);
            OPTS.BC(1) = 0;                 % Set electric wall at left
        otherwise
            if(VMODE>0), fprintf('unknown (''%s''). ERROR.\n', OPTS.hsymm); end
            error('[MP] sisolver3d - unknown symmetry condition specified');
    end
end


% For bend modes:
if (nargin > 5) && isfield(OPTS,'radius');
    R = OPTS.radius;  N.x = N.x + R;      % Set left edge of domain to be at radius
    if(VMODE>0), fprintf([mfilename ': Using radius of %g um at LEFT EDGE of the computational domain.\n'], R); end
%    mu_guess = mu_guess * R;              % Reset eigenvalue
    OPTS.coordmode = 'R';
%    PMLwidth = [0 0.5 0 0]; PMLsigma = [0.1 0.1]; % um, um^-1
end

% Call modesolver engine with the prepared parameters
% if(isfield(OPTS,'enginever'))
if(VMODE>0), disp(['Calling solver engine ' OPTS.enginever '.']); end
[beta, F, V] = feval(OPTS.enginever, N, k0, mu_guess, OPTS, NMODES_CALC, OPTS.PMLwidth, OPTS.PMLsigma);   % Call custom solver
% else
%     if(VMODE>0) disp('Calling default engine (m2wcyl).'); end
%     [beta, F, V] = m2wcyl(N, k0, mu_guess, OPTS, NMODES_CALC, OPTS.PMLwidth, OPTS.PMLsigma);                  % Call default solver
% end

% Compute auxiliary loss value outputs from output eigenvalues
if(strcmp(OPTS.coordmode,'R'))      % For circularly bent structures (cylindrical)
    alpha = imag(beta)/R;           % beta here is really gamma, were mode has exp(i*gamma*phi) field dependence
    F.dB90 = 20*alpha*(pi*R/2)*log10(exp(1));       % This does not depend on choice of length units of inputs
%     F.dB90 = 20*imag(beta)*(pi/2)*log10(exp(1));    % This does not depend on choice of length units of inputs
else                                % For straight structures (Cartesian)
    alpha = imag(beta);
end
F.dBcm = 0.2/log(10) * alpha;       % Loss in dB/cm if input lengths are in meters
                                    % LossdB = -10*log10( P(z=1 unit)/P(z=0) ) = -10*log10( exp(-2*alpha*1 unit) ) = 20*alpha*1 unit/ln(10)

% Store various useful data in output data structures
F.beta = beta;
%warning off MATLAB:divideByZero;
F.LossQ = real(beta)./imag(beta)/2; % Rough LossQ ~ beta / 2 alpha; true correct value should be LossQ = k0 ng / 2 alpha (i.e. with group index, ng)
%warning on MATLAB:divideByZero;

% Remove undesired modes
F.Ex = F.Ex(:,:,1:NMODES_KEEP); F.Ey = F.Ey(:,:,1:NMODES_KEEP); F.Ez = F.Ez(:,:,1:NMODES_KEEP);
F.Hx = F.Hx(:,:,1:NMODES_KEEP); F.Hy = F.Hy(:,:,1:NMODES_KEEP); F.Hz = F.Hz(:,:,1:NMODES_KEEP);
F.dx = dx; F.dy = dy;
%Fb.Pr = Pr; Fb.Pz = Pz; Fb.Rmtx = Rmtx;
end