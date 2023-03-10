%---------.---------.---------.---------.---------.---------.---------.--------|
% (under construction)
% June 18, 2002 Milos Popovic
% Cylindrical-PML 2D mode solver for arbitrary index distributions.
%
% function [mu, F, {V, D, Pr, Pz, Rmtx}] = m2wcyl(S, nu_index, mu_guess, OPTS, nmodes, dPML, sigmaPML)
%
% ** WARNING **    Modal fields returned [Ex Ey Ez Hx Hy Hz] are defined
%                  along a left-handed Cartesian coordinate system, in
%                  order to stay consistent with rho-phi-z cylindrical
%                  coordinates. I know this is stupid but using [Ex Ez
%                  Ey..] would be just as confusing. Thus if you are using
%                  Cartesian coords, on a right-handed (standard) coord
%                  system, the field values are really [Ex -Ey Ez Hx -Hy
%                  Hz].
%
% OUTPUT:
% mu             = ('w'-mode) energy-mode complex resonant frequencies, w/c (highest Q first, then lowest frequency first)
%                  ('b'-mode) power-mode complex propagation constants, gamma =def= beta*R (lowest loss first, then highest eff. index first)
%                  NOTE: for adjoint computations, modes are sorted by largest loss first
% F              = data structure containing fields of computed modes
% V,D            = exact eigenvalue (D)/eigenvector (V) matrices returned by eigs without post-processing
% Pr, Pz, Rmtx   = PML complex coordinate stretching factors Pr, Pz and complex radial coordinate Rmtx
%
% INPUT:           ("pixel" grid is (M-1)x(N-1))
% S              = structure: S.n, S.x(:), S.y(:) = (2M-1 x 2N-1, 1x2M-1, 1x2N-1)
% nu_index       = index of eigenvalue eqn:
%                   # lambda's around in w-mode, frequency w/c in b-mode
% mu_guess       = eigenvalue guess: w guess in w-mode, beta*R (gamma) guess in b-mode
%    NOTE: mu_guess is converted to w/c ('w'-mode) or gamma (beta*R = beta in cartesian coords, 'b'-mode) inside the function!
% OPTS.coordmode = cylind'R'ical or 'C'artesian
% OPTS.eigmode   = 'w' or 'b' eigenvalue
% OPTS.fieldmode = 'V'ector or se'MX'ivector (or 'MY' for y-polarization)
%                    'V'ector mode solves complete field problem; 'MX'
%                    solves only for Ex (horizontal) field while 'MY'
%                    solves only for Ey (vertical) field, disregarding
%                    Ex-Ey coupling (reduces problem by 2x in size).
% OPTS.BC        = [OPTIONAL] default = [0 0 0 0]; [left right bottom top] boundary conditions,
%                    PEC = 0; PMCh = 1; PBC = 2; PMC = 3; PECh = 4.
%                    PEC  = perfect electric conductor, defined at pixel edge
%                    PMC  = perfect magnetic conductor, defined at pixel edge
%                    PECh = perfect electric conductor, defined at pixel center
%                    PMCh = perfect magnetic conductor, defined at pixel center
%                    PBC  = periodic boundary condition
% OPTS.adjoint   = [OPTIONAL] default = 0; 0 = normal ([Ex Ey]), 1 = adjoint ([Hy -Hx])
%                    Return solution for adjoint operator.
% OPTS.eigsfcn   = [OPTIONAL] default = 'eigs'; can use to set custom sparse eigenvalue solver
%                    (Generally internal use) Sets which Matlab function is
%                    used to solve the sparse eigenvalue problem.
% OPTS.operver   = [OPTIONAL] default = 'm2dpmloperR2'; can also set to old version m2dpmloper
% nmodes         = number of modes to find/compute
% dPML           = [4-vector, Units: length] PML thickness: [left right down up]
%                    The PML (perfectly matched layer) is a medium which reaches
%                    into the computational domain from the domain boundaries and
%                    retains the underlying refractive index but has a quadratically
%                    tapered attenuation from zero to maximum at the boundary. e.g. [0 0.5 0.5 0] (microns)
% sigmaPML       = [2-vector, Units: 1/length] Maximum PML conductivity: [horizontal vertical]
%                    Maximum conductivity of the PML medium (i.e. the value
%                    of the conductivity at the boundaries). e.g. [0.2 0.2] (1/microns)
%                    The simulations are insensitive to this parameter (can
%                    be +/- a couple orders of magnitude and PML will work).


% Known Issues: - Omega_tilde, Gamma_tilde/Gamma_hat offsets not yet
% -------------   implemented for proper phase offsets in discrete time/z (FDTD).
%               - Must fix H-field (maybe also E-field) computation at the
%                 end to work for PMC boundary conditions!
%               - [Nov272003] Need to fix 'MY' semivector mode to be able to keep only
%                 one-pixel-tall domain for TE slab calculations (TM already works).
% Other todo:   - Create a few-page Latex PDF manual giving the math, PML
% -----------     definitions, interface and simple examples, and...
%                   - Explain V,MX,MY field/polarization modes.
%                   - Explain boundary condition types.
%
% Updates:
% --------
% Jun 18, 2002 - First version.
% Sep 22, 2002 - Completing ability to compute for beta-eigenvalues.
%              - NOTE: Ez is still *always* calculated assuming cartesian
%                coordinates! This needs to be fixed.
%              - Lparam (nu_index) definition changed from w^2 to k0 for 'b'-mode operation.
% Oct  1, 2002 - Lparam renamed to nu_index; eigguess renamed to mu_guess.
%              - Set VMODE = 1 to display textual debug info; VMODE = 2 plots figures also (all mode field plots..)
% Oct  3, 2002 - Changed EOPTS options structure for "eigs" to OPTS so that options for eigs can be passed IN to m2wcyl, such as
%                a starting eigenvector. Raw V and D (eigenvalues and eigenvectors) are also an available output).
% Nov 20, 2002 - Proper calculation of third E-field (Ez, i.e. Ephi) added; H-field calculation?
% Nov 26, 2002 - Proper calculation of H-field (for radial or Cartesian coords) added.
% Dec 11, 2002 - Added optional parameter OPTS.adjoint = 0 (default) or 1 (adjoint),
%                allowing one to compute the modes of the adjoint system.
%                *NOTE* that for modes of the adjoint system, [Ex Ey] -> [Hy -Hx], ie.
%                the output fields are different ones.
% Feb 23, 2003 - MOPTS changed to match new version of m2dpmloper, to pass
%                fieldmode/eigmode instead of fmode/emode.
% Jun  7, 2003 - Added Pr, Pz, Rmtx as optional output variables (for leaky mode
%                overlap integrals without adjoint).
%              - FIXED calculation of Ez from Ex and Ey to be properly
%                computed using Pr, Pz in the PML; also now works in
%                'w'-mode! The E-field should now be exactly correct. Still
%                need to fix H calculation, at least for 'w'-mode (check).
% Jun  9, 2003 - FIXED to pass OPTS not MOPTS to m2dpmloper so that BC's
%                get passed. This was fixed in the /workFeb16/ directory version of
%                m2wcyl, so that version and this one must be made consistent ASAP!!
% Jul 23, 2003 - ADDED ability to pass in a user-set "eigs" function in the
%                OPTS.eigsfcn variable.
% Jul 28, 2003 - *** MAJOR UPDATE & SYNCHRONIZATION BETWEEN MULTIPLE
%                VERSIONS OF m2wcyl.m ***
%              --Inserted updates of other m2wcyl versions:--
%                **Jun  9, 2003 - quick fix: can now pass on BCs (done in other versions of
%                m2wcyl better).
%                Jun 25, 2003 - FIXED calculation of Hx, Hy, Hz fields for beta-mode for
%                cylindrical, leaky modes. (still need to do for
%                omega-modes).
%              -----------------------------------------------
%              - Renamed w to mu
%              - Cleaned up commented out lines
%              - **Changed input argument mu_guess from w-units to have
%                k0-units (w/c)
%              - **Changed wI to be <0 as it should be for loss!
%              - Verified and cleaned up output field calculations
%              - Fixed Omega_tilde to have negative sign (this isn't used yet though)
%              - Corrected H-field computation in 'w'-mode (the guess
%                frequency was being used!! must use computed complex frequency).
% Aug  1, 2003 - Corrected H-field computation for PMC boundary conditions.
%                Still have to check whether this works for cylindrical modes!
% Aug 11, 2003 - Fixed H-field computation to include PML factors (Pr, Pz).
%                It should now be exactly right everywhere (domain + PML).
% Oct 23, 2003 - Added option OPTS.sigma = [] to search for closest
%                eigenvalue (DEFAULT), or 'SR','SI','LR', ...
% Oct 27, 2003 - Added default setting of OPTS.BC to [0 0 0 0] in case it
%                is not set at input. This was done just for compatibility
%                with older script files.
% Nov 21, 2003 - Updates: Adding field capture for PMCh BCs (different size)
% Mar 21, 2005 - Minor edit to set OPTS.tol = eps *only* if user doesn't supply it
%              - Modified to permit multiple guess values in mu_guess - the
%                operator is normalized by mu_guess(1), as is the vector
%                mu_guess.
% Aug 30, 2006 - Minor fix to add abs() to allow imagesc of complex index
%                distributions when VMODE is on.
% Nov     2006 - IN PROGRESS: Adding Star-P parallelism capability
% Nov  3, 2006 - IN PROGRESS: Modifying to support new recoded m2dpmloperR2
%              - Added OPTS.operver to allow selection between old (m2dpmloper)
%                and new (m2dpmloperR2) operator function.
% Nov  4, 2006 - Added PECh definition
%              - Created derivative function for doing gradients and curls
%                of all fields with boundary conditions properly considered
% Dec  9, 2010 - Commented out "pack" function because new versions of Matlab
%                don't support it inside functions.
% Sep  8, 2011 - Help text at top of file, OPTS.operfcn corrected to be
%                OPTS.operver (the right name used in fcn below).
% Apr 30, 2012 - Updated help text to clarify variables and their meaning
%              - Changed default value of OPTS.operver to 'm2dpmloperR2'
%              - Went through code to make it conform to Matlab R2011b standard
%

function [mu, F, V, D, Pr, Pz, Rmtx] = m2wcylR3(S, nu_index, mu_guess, OPTS, nmodes, dPML, sigmaPML)
PEC = 0; PMCh = 1; PBC = 2; PMC = 3; PECh = 4;          % [MP-21nov2003] Define boundary condition identifiers
c = 299792458; u0 = 4e-7*pi; %e0 = 1/c^2/u0;            % Physical Constants
global VMODE;  if isempty(VMODE), VMODE = 0; end        % [VMODE] Default value for debug flag.
if ~isfield(OPTS,'adjoint'), OPTS.adjoint = 0; end      % Default "adjoint-computation" flag is FALSE (0).
if isempty(OPTS.adjoint), OPTS.adjoint = 0; end
if ~isfield(OPTS,'eigsfcn'), OPTS.eigsfcn = 'eigs'; end % Default sparse-matrix eigenvalue solver function.
if isempty(OPTS.eigsfcn), OPTS.eigsfcn = 'eigs'; end
if ~isfield(OPTS,'sigma'), OPTS.sigma = []; end         % Eigenvalue choice: shift-invert = [], or 'LM','SI',etc
if isempty(OPTS.sigma), OPTS.sigma = mu_guess/mu_guess(1); end    % Guess not needed for direct modes.mu_guess = 1; (for multiple eigenvalues)
if ~isfield(OPTS,'BC'), OPTS.BC = [0 0 0 0]; end        % Default boundary conditions
if isempty(OPTS.BC), OPTS.BC = [0 0 0 0]; end
BCl = OPTS.BC(1); BCr = OPTS.BC(2); BCd = OPTS.BC(3); BCu = OPTS.BC(4); % Macros [MP-21nov2003]
if ~isfield(OPTS,'operver'), OPTS.operver = 'm2dpmloperR2'; end % Default operator function
if isempty(OPTS.operver), OPTS.operver = 'm2dpmloperR2'; end

t0 = clock;                                             % Begin timing the code (for debug)

R = S.x; Z = S.y; nn = S.n;                             % [MPREM] Import parameters for now...
L = nu_index; % l0 = mu_guess;                          % [MPREM] [MP] removed Sep 22, 2002
sigmaRMAX = sigmaPML(1); sigmaZMAX = sigmaPML(2);

if VMODE > 1                                            % [VMODE] Plot refr. index distribution
%   figure; imagesc(R, Z, (nn.' - min(nn(:))) * 1000); axis image;
   figure; imagesc(R, Z, abs(nn.' - min(nn(:))) * 1000); axis image;    % [20060830-MP] Include abs() for complex indices
   if (OPTS.coordmode == 'C')
       xlabel('x'); ylabel('y');
   elseif (OPTS.coordmode == 'R')
       xlabel('\rho'); ylabel('z');
   end
   title('2D Index Profile'); filestampplot(gcf, mfilename, 'I'); % Label w/ filename
end

NR = (length(R)+1)/2; NZ = (length(Z)+1)/2;
if OPTS.eigmode == 'w'      % For 'w'-eigenvalue case...
disp('Warning: new version of modesolver uses mu_guess = k0 = w/c as frequency guess value, NOT w as before\n.');
%    mu_guess = mu_guess/c;  % [MP] need to change the guess value later to k0! (now it is w0)
    w0 = c*mu_guess(1);     % Modified in case a vector of guess values is supplied (possible for some eigensolvers, e.g. jdqr.m)
else                        % For 'beta'-eigenvalue case...
    w0 = c*nu_index;        % [MP] Sep 20, 2002
end


% Generate PML matrices
[Pr, Pz, Rmtx] = m2dpmlmatx(w0, nn.^2, R, Z, [sigmaRMAX sigmaZMAX], dPML);  % w0 is exact in 'b'-mode or the guess value in 'w'-mode

% Generate "vector-Helmholtz" operator matrix for 'w'- or 'b'-eigenvalue problem
dR = R(3) - R(1); dZ = Z(3) - Z(1);                     % Grid spacing in R, Z is dR/2, dZ/2!
if (OPTS.coordmode == 'C'), Rmtx = ones(size(Rmtx)); end % [MP] Oct 01, 2002 - For Cartesian mode, set Rmtx=all 1's.
%H = m2dpmloper(nn.^2, Pr, Pz, Rmtx, dR, dZ, L^2, OPTS); % Works for 'w'- and 'b'-mode operation
H = feval(OPTS.operver, nn.^2, Pr, Pz, Rmtx, dR, dZ, L^2, OPTS); % Works for 'w'- and 'b'-mode operation; modified call (feval) allows arbitrary operator fcn to be called
H = H / mu_guess(1)^2;                                     % Normalize matrix to have modal (k/k0)^2 as eigenvalue in 'w'-mode
%pack;                                                   % [MPREM] MEMORY Pack.

t1 = clock;

% Field distribution boundary offset indices depending on boundary conditions
%oxL = 1*(BCl == PMCh); oxR = 1*(BCr == PMCh); oxD = 0; oxU = 0; % [MP-21nov2003] Index offsets for Ex field
%oyL = 0; oyR = 0; oyD = 1*(BCd == PMCh); oyU = 1*(BCu == PMCh); % [MP-21nov2003] Index offsets for Ey field
%MxSIZE = (NR-1-(oxL+oxR))*(NZ-2-(oxD+oxU)); MySIZE = (NR-2-(oyL+oyR))*(NZ-1-(oyD+oyU)); % Size of Ex and Ey parts of the solution vector
oxL = 0; oxR = 0; oxD = 1*(BCd == PEC || BCd == PECh); oxU = 1*(BCu == PEC || BCu == PECh || BCu == PBC); % [2006-11-04-MP] Index offsets for Ex field matrix size reduction for some BCs
oyL = 1*(BCl == PEC || BCl == PECh); oyR = 1*(BCr == PEC || BCr == PECh || BCr == PBC); oyD = 0; oyU = 0; % [2006-11-04-MP] Index offsets for Ey field matrix size reduction for some BCs
MxSIZE = (NR-1-(oxL+oxR))*(NZ-(oxD+oxU)); MySIZE = (NR-(oyL+oyR))*(NZ-1-(oyD+oyU)); % Size of Ex and Ey parts of the solution vector
if OPTS.fieldmode(1) == 'M'                             % Reduce matrix to single polarization!
    if OPTS.fieldmode(2) == 'X'                         % Semi-vectorial operator for Ex modes
        H = H(1:MxSIZE,1:MxSIZE);
    elseif OPTS.fieldmode(2) == 'Y'                     % Semi-vectorial operator for Ey modes
        H = H(MxSIZE+(1:MySIZE),MxSIZE+(1:MySIZE));
    else
        error([mfilename ': For semi-vectorial calculations, solvertype must '...
               'be ''MX'' or ''MY''.']);
    end
end

if (OPTS.adjoint), H = H'; end;                         % For adjoint system computation, take adjoint of operator.
asgn = (-1)^OPTS.adjoint;                               % Sign flag used for sorting eigenvalues of the normal (1) and adjoint (-1) system.

% Allow flag to permit StarP parallel mode solving
if(isfield(OPTS,'PARALLEL') && OPTS.PARALLEL),  H = matlab2pp(H);  end


% Solve eigenvalue problem by Arnoldi algorithm in shift-invert mode
if(~isfield(OPTS,'tol')),  OPTS.tol = eps;  end;  %1e-10;   % [MP] Tolerance set *only* if not passed in by user in OPTS structure
OPTS.disp = VMODE;
%[V,D] = eigs(H,nmodes,1,OPTS);                         % normal Matlab solver
%[V,D] = feval(OPTS.eigsfcn, H, nmodes, 1, OPTS);        % general Matlab eigs solver (allows user-set eigs function)
[V,D] = feval(OPTS.eigsfcn, H, nmodes, OPTS.sigma, OPTS);        % general Matlab eigs solver (allows user-set eigs function)
%[V,D] = feval(OPTS.eigsfcn, H, nmodes, 'LR', OPTS);        % general Matlab eigs solver (allows user-set eigs function)
%pack;                                                   % [MPREM]
mu = sqrt(diag(D)) * mu_guess(1);                          % mu eigenvalue = w/c ('w'-mode) or gamma ('b'-mode)

% Sort solved modes by lowest imaginary part of propagation constant!
if (OPTS.eigmode == 'w')    % For 'w'-eigenvalue case...
%    X = sortrows([-imag(w)*asgn -real(w) V.']);   % [MPREM] -ve imaginary part since exp(-i*w*t): w =def= wR - i wI  (1/tau =def= wI)
    X  = sortrows([real(mu) +imag(mu)*asgn V.'],'descend');      % -ve imaginary part since exp(-i*w*t): w =def= wR - i wI  (1/tau =def= wI)
            % Sort second by LOWEST frequency first.
%    mu = +X(:,2) + i*X(:,1)*asgn;  % [MPREM] swap columns wrt top part of this "if"-statement.
    mu = X(:,1) + 1i*X(:,2)*asgn;  % keeping wI < 0, swap columns wrt top part of this "if"-statement.
else                        % ...and for beta-eigenvalue case.
    % If betas have imaginary parts, sort by lowest loss (highest in adjoint case)...
    % ... otherwise, sort by highest modal index (best confined mode first)
    X = sortrows([real(mu) +imag(mu)*asgn V.'],'descend');       % +ve imaginary part since exp(+i*beta*z): beta =def= bR + i bI  (alpha =def= bI)
    mu = X(:,1) + 1i*X(:,2)*asgn;  % swap columns wrt top part of this "if"-statement.
end
V = X(:,3:size(X,2)).';                                 % [MP030728] Resort V,D raw matrices to match above sorting.
D = diag( (mu/mu_guess).^2 );


t2 = clock;                                             % End timing the eigenvalue solver (ARPACK)
if VMODE ~= 0                                           % Print elapsed time
    fprintf(['Timing of the code:\n' ...
             '---------------------------------------------\n']);
    fprintf('"Vector-Helmholtz" operator setup: %5.2f sec\n',etime(t1,t0));
    fprintf('ARPACK (double) eigenvalue solver: %5.2f sec\n',etime(t2,t1));
    fprintf('---------------------------------------------\n');
    fprintf('TOTAL:                             %5.2f sec\n',etime(t2,t0));
    if VMODE > 1
        figure; plot( real(sqrt(diag(D)))   ); title('Normalized eigenvalue real parts');
        xlabel('Mode number'); ylabel('Re(\mu)');       % Plot real part of mu (=neff*k0 in 'b'-mode case)
    end
end


% OUTPUT MODE FIELDS
if nargout > 1                                          % Extracting mode patterns from eigenvector
    NX = NR; NY = NZ;                                   % [MPREM] [MP] temporary..
%    F.Rr = R(2:2:2*NR-2); F.Zr = Z(3:2:2*NZ-3);
%    F.Rz = R(3:2:2*NR-3); F.Zz = Z(2:2:2*NZ-2);
    F.Rr = R((2+2*oxL):2:(2*NR-2-2*oxR)); F.Zr = Z((1+2*oxD):2:(2*NZ-1-2*oxU));     % Set up 
    F.Rz = R((1+2*oyL):2:(2*NR-1-2*oyR)); F.Zz = Z((2+2*oyD):2:(2*NZ-2-2*oyU));

    er = nn.^2;                                         % [MPREM]
    % Preallocate field matrices
    Ex = zeros(NX-1-(oxL+oxR),NY-(oxD+oxU),size(V,2));
    Ey = zeros(NX-(oyL+oyR),NY-1-(oyD+oyU),size(V,2));
    betaR = zeros(1,size(V,2));
    omega = zeros(1,size(V,2));
    Gammah= zeros(1,size(V,2)); Gammat= zeros(1,size(V,2));
    Omegat= zeros(1,size(V,2));
    for k = 1:size(V,2)                                 % Extract fields from eigenvectors ...
%         if(OPTS.fieldmode == 'V')                       % ... based on fieldmode in solver call.
%             % [MP-21nov2003] Extract sizes according to boundary condition type!
%             Ex(:,:,k) = [zeros(oxL,NY-2); zeros(NX-1-(oxL+oxR),oxD), ...
%                          reshape(V(1:MxSIZE,k),NX-1-(oxL+oxR),NY-2-(oxD+oxU)), ...
%                          zeros(NX-1-(oxL+oxR),oxU); zeros(oxR,NY-2)];
%             Ey(:,:,k) = [zeros(oyL,NY-1); zeros(NX-2-(oyL+oyR),oyD), ...
%                          reshape(V(MxSIZE+1:MxSIZE+MySIZE,k),NX-2-(oyL+oyR),NY-1-(oyD+oyU)), ...
%                          zeros(NX-2-(oyL+oyR),oyU); zeros(oyR,NY-1)];
%         elseif(OPTS.fieldmode == 'MX')
%             Ex(:,:,k) = [zeros(oxL,NY-2); zeros(NX-1-(oxL+oxR),oxD), ...
%                          reshape(V(1:MxSIZE,k),NX-1-(oxL+oxR),NY-2-(oxD+oxU)), ...
%                          zeros(NX-1-(oxL+oxR),oxU); zeros(oxR,NY-2)];
%             Ey(:,:,k) = zeros(NX-2,NY-1);
%         elseif(OPTS.fieldmode == 'MY')
%             Ex(:,:,k) = zeros(NX-1,NY-2);
%             Ey(:,:,k) = [zeros(oyL,NY-1); zeros(NX-2-(oyL+oyR),oyD), ...
%                          reshape(V(1:MySIZE,k),NX-2-(oyL+oyR),NY-1-(oyD+oyU)), ...
%                          zeros(NX-2-(oyL+oyR),oyU); zeros(oyR,NY-1)];
%         else
%             error('m2dchew: solvertype invalid; must be ''V'', ''MX'' or ''MY''.');
%         end
        switch OPTS.fieldmode
            case {'V','MX'}
                % [2006-11-04-MP] Extract field matrices with sizes determined by boundary condition type
                Ex(:,:,k) = reshape(V(1:MxSIZE,k),NX-1-(oxL+oxR),NY-(oxD+oxU));
            case 'MY'
                Ex(:,:,k) = zeros(NX-1-(oxL+oxR),NY-(oxD+oxU));
            otherwise
                error([mfilename ': solver fieldmode invalid; must be ''V'', ''MX'' or ''MY''.']);
        end
        switch OPTS.fieldmode
            case {'V'}
                Ey(:,:,k) = reshape(V(MxSIZE+(1:MySIZE),k),NX-(oyL+oyR),NY-1-(oyD+oyU));
            case {'MY'}
                Ey(:,:,k) = reshape(V(1:MySIZE,k),NX-(oyL+oyR),NY-1-(oyD+oyU));
            case 'MX'
                Ey(:,:,k) = zeros(NX-(oyL+oyR),NY-1-(oyD+oyU));
            otherwise
                error([mfilename ': solver fieldmode invalid; must be ''V'', ''MX'' or ''MY''.']);
        end

        % Set propagation constant along phi/"z" direction
        if (OPTS.eigmode == 'b')
            betaR(k) = mu(k); omega(k) = w0;
        else
            betaR(k) = nu_index; omega(k) = c*mu(k);    % [MP030728] Fixed field derivation for 'w'-mode
        end
        % Derive Ez from Gauss' Law
        % (here imaginary unit 'i' uses the exp(-i w t), exp(+i beta z)
        % convention (physics convention). When EE j is used, it is
        % explicitly written as j, NEVER as i below.
%        erxx = er(2:2:2*NX-2,3:2:2*NY-3); eryy = er(3:2:2*NX-3,2:2:2*NY-2);
%        erzz = er(3:2:2*NX-3,3:2:2*NY-3);
%        rxx = Rmtx(2:2:2*NX-2,3:2:2*NY-3); ryx = Rmtx(3:2:2*NX-3,3:2:2*NY-3);
        jxx = (1+2*oxD):2:(2*NY-1-2*oxU); iyy = (1+2*oyL):2:(2*NX-1-2*oyR);
        erxx = er(2:2:2*NX-2,jxx); eryy = er(iyy,2:2:2*NY-2);     % [2006-11-04] Corrected for new BCs
        erzz = er(iyy,jxx); rxx = Rmtx(2:2:2*NX-2,jxx); ryx = Rmtx(iyy,jxx);

%        Ax = rxx .* (erxx).*[Ex(:,:,k)]; Ax = Pr(3:2:2*NX-3,3:2:2*NY-3) .* (Ax(2:NX-1,:) - Ax(1:NX-2,:))/dR;  % d/dx exx Ex
%        Ay = (eryy).*[Ey(:,:,k)]; Ay = Pz(3:2:2*NX-3,3:2:2*NY-3) .* ryx .* (Ay(:,2:NY-1) - Ay(:,1:NY-2))/dZ;  % d/dy eyy Ey
        Ax = rxx .* erxx .* Ex(:,:,k);
%         if(BCl == PMC)  Ax = [-Ax(1,:); Ax];                               % Add left field point for gradient only for PMC, PMCh or PBC
%             elseif(BCl == PMCh) Ax = [zeros(1,size(Ax,2)); Ax];
%             elseif(BCl == PBC)  Ax = Ax([end 1:end],:); end
%         if(BCr == PMC)  Ax = [Ax; -Ax(end,:)];                             % Add right field point for gradient only for PMC and PMCh
%             elseif(BCr == PMCh) Ax = [Ax; zeros(1,size(Ax,2))]; end
%         Ax = Pr(iyy,jxx) .* diff(Ax,1,1)/dR;  % d/dx exx Ex
         Ax = Pr(iyy,jxx) .* diffxy(Ax,'Ex',[BCl BCr BCd BCu],0)/dR;  % d/dx exx Ex

        Ay = eryy .* Ey(:,:,k);
%         if(BCd == PMC)  Ay = [-Ay(:,1), Ay];                               % Add bottom field point for gradient only for PMC, PMCh or PBC
%             elseif(BCd == PMCh) Ay = [zeros(size(Ay,1),1), Ay];
%             elseif(BCd == PBC)  Ay = Ay(:,[end 1:end]); end
%         if(BCu == PMC)  Ay = [Ay, -Ay(:,end)];                             % Add top field point for gradient only for PMC and PMCh
%             elseif(BCu == PMCh) Ay = [Ay, zeros(size(Ay,1),1)]; end
%         Ay = Pz(iyy,jxx) .* ryx .* diff(Ay,1,2)/dZ;  % d/dy eyy Ey
         Ay = Pz(iyy,jxx) .* ryx .* diffxy(Ay,'Ey',[BCl BCr BCd BCu],1)/dZ;  % d/dy eyy Ey
        
        deltaphi = 0; Gammah(k) = betaR(k) * exp(-1i*deltaphi);  % Define delta-phi and Gamma-hat-prime; [MP030728] Verified eqn
%        warning off MATLAB:divideByZero                         % [MP] Added July 26, 2003
        if( Gammah(k) == 0 ), disp('Warning: Ez improper as betaR = 0; Ez should be =def= 0 here.'); end % [MP030728]
        Ez(:,:,k) = 1i./(Gammah(k)*erzz) .* (Ax + Ay);           % [MP020728] Verified eqn
%        warning on  MATLAB:divideByZero
        clear erxx eryy erzz;                                   % [MP030728] Added for memory conservation


        % Derive H fields from Faraday's Law  ([MP030728] eqns verified again)
        if(OPTS.adjoint), disp('Warning: computation of other 4 field components not supported in "adjoint" mode.'); end;  % [MP] Jul 15, 2004
        deltat = 0; Omegat(k) = omega(k) * exp(-1i*deltat);     % Defining delta_t and Omega_tilde, CHANGED SIGN [MP030728]
        C = -1i/(Omegat(k) * u0);
        Gammat(k) = betaR(k) * exp(+1i*deltaphi);               % Define delta-phi and Gamma-tilde, [MP030728]
%        ryy = Rmtx(3:2:2*NX-3,2:2:2*NY-2); rzz = Rmtx(3:2:2*NX-3,3:2:2*NY-3);
        ryy = Rmtx(iyy,2:2:2*NY-2); rzz = Rmtx(iyy,jxx);                   % [2006-11-04-MP]
% %        Hx(:,:,k) = C * (i*Gammat(k).*Ey(:,:,k)./ryy - ...
% %                    ([Ez(:,:,k) zeros(NX-2,1)] - [zeros(NX-2,1) Ez(:,:,k)])/dZ );
% %        Hz(:,:,k) = C * (([Ex(:,:,k) zeros(NX-1,1)] - [zeros(NX-1,1) Ex(:,:,k)])/dZ - ...
% %                    ([Ey(:,:,k); zeros(1,NY-1)] - [zeros(1,NY-1); Ey(:,:,k)])/dR);
% %        Hy(:,:,k) = C * (1./rxx .* ...
% %                    ([rzz .* Ez(:,:,k); zeros(1,NY-2)] - [zeros(1,NY-2); rzz .* Ez(:,:,k)])/dR - ...
% %                   i * Gammat(k) * Ex(:,:,k)./rxx);
% %       Taking care of BCs... Recall: PEC = 0; PMC = 1; % IF MORE BCS ARE
% %       TO BE ADDED, MUST CHANGE LINES BELOW TO COMPUTE SEPARATELY BASED ON
% %       WHICH BC IS ACTIVE...
% %        Hx(:,:,k) = C * (i*Gammat(k).*Ey(:,:,k)./ryy - ...
% %                    ([Ez(:,:,k) Ez(:,NY-2,k)*OPTS.BC(4)] - [Ez(:,1,k)*OPTS.BC(3) Ez(:,:,k)])/dZ );
% %        Hz(:,:,k) = C * (([Ex(:,:,k) Ex(:,NY-2,k)*OPTS.BC(4)] - [Ex(:,1,k)*OPTS.BC(3) Ex(:,:,k)])/dZ - ...
% %                    ([Ey(:,:,k); Ey(NX-2,:,k)*OPTS.BC(2)] - [Ey(1,:,k)*OPTS.BC(1); Ey(:,:,k)])/dR);
% %        Hy(:,:,k) = C * (1./rxx .* ...
% %                    ([rzz .* Ez(:,:,k); rzz(NX-2,:) .* Ez(NX-2,:,k) * OPTS.BC(2)] - ...
% %                     [rzz(1,:) .* Ez(1,:,k) * OPTS.BC(1); rzz .* Ez(:,:,k)])/dR - ...
% %                    i * Gammat(k) * Ex(:,:,k)./rxx);
%         Hx(:,:,k) = C * (i*Gammat(k).*Ey(:,:,k)./ryy - ...
%                     Pz(3:2:2*NX-3,2:2:2*NY-2) .* ([Ez(:,:,k) Ez(:,NY-2,k)*OPTS.BC(4)] - [Ez(:,1,k)*OPTS.BC(3) Ez(:,:,k)])/dZ );
%         Hz(:,:,k) = C * (Pz(2:2:2*NX-2,2:2:2*NY-2) .* ([Ex(:,:,k) Ex(:,NY-2,k)*OPTS.BC(4)] - [Ex(:,1,k)*OPTS.BC(3) Ex(:,:,k)])/dZ - ...
%                     Pr(2:2:2*NX-2,2:2:2*NY-2) .* ([Ey(:,:,k); Ey(NX-2,:,k)*OPTS.BC(2)] - [Ey(1,:,k)*OPTS.BC(1); Ey(:,:,k)])/dR);
%         Hy(:,:,k) = C * (1./rxx .* ...
%                     Pr(2:2:2*NX-2,3:2:2*NY-3) .* ([rzz .* Ez(:,:,k); rzz(NX-2,:) .* Ez(NX-2,:,k) * OPTS.BC(2)] - ...
%                      [rzz(1,:) .* Ez(1,:,k) * OPTS.BC(1); rzz .* Ez(:,:,k)])/dR - ...
%                     i * Gammat(k) * Ex(:,:,k)./rxx);

        Hx(:,:,k) = C * (1i*Gammat(k).*Ey(:,:,k)./ryy - ...
                    Pz(iyy,2:2:2*NY-2) .* diffxy(Ez(:,:,k),'Ez',[BCl BCr BCd BCu],1)/dZ);
        Hz(:,:,k) = C * (Pz(2:2:2*NX-2,2:2:2*NY-2) .* diffxy(Ex(:,:,k),'Ex',[BCl BCr BCd BCu],1)/dZ - ...
                    Pr(2:2:2*NX-2,2:2:2*NY-2) .* diffxy(Ey(:,:,k),'Ey',[BCl BCr BCd BCu],0)/dR);
        Hy(:,:,k) = C * (1./rxx .* ...
                    Pr(2:2:2*NX-2,jxx) .* diffxy(rzz .* Ez(:,:,k),'Ez',[BCl BCr BCd BCu],0)/dR - ...
                    1i * Gammat(k) * Ex(:,:,k)./rxx);


        if( VMODE > 1 )                                 % For guided modes, Ex/Ey=real, Ez=imag (90? out of phase!).
            V1=Ex(:,:,k); i1=find(abs(V1) == max(abs(V1(:)))); A1 = V1(i1(1));
            figure; imagesc(F.Rr, F.Zr, real(Ex(:,:,k).'/A1)); axis image; % E^X
            title(strcat('E^X (Mode ',num2str(k),')')); xlabel('x'); ylabel('y');

            V1=Ey(:,:,k); i1=find(abs(V1) == max(abs(V1(:)))); A1 = V1(i1(1));
            figure; imagesc(F.Rz, F.Zz, real(Ey(:,:,k).'/A1)); axis image; % E^Y
            title(strcat('E^Y (Mode ',num2str(k),')')); xlabel('x'); ylabel('y');

            V1=Ez(:,:,k); i1=find(abs(V1) == max(abs(V1(:)))); A1 = V1(i1(1));
            figure;	imagesc(F.Rz, F.Zz, real(Ez(:,:,k).'/A1)); axis image; % E^Z
            title(strcat('E^Z (Mode ',num2str(k),')')); xlabel('x'); ylabel('y');
            clear V1 i1;                                % [MP030728] Added for memory conservation
        end
    end
end

% Display elapsed time
if( VMODE ~= 0 ),  fprintf('Elapsed time: %6.2f sec\n', etime(clock,t0));  end

% mu already defined; putting fields into data structure F
F.Ex = Ex; F.Ey = Ey; F.Ez = Ez; F.Hx = Hx; F.Hy = Hy; F.Hz = Hz;
F.BC = [BCl BCr BCd BCu];                                                  % [2006-11-04] Adding boundary conditions to fields data structure
F.ver = 'R2';   % Stamp solver version into fields structure to allow other codes to recognize data format
%---------.---------.---------.---------.---------.---------.---------.--------|


% Field differentiation helper functions which know boundary conditions
% and mirror the fields appropriately around the boundaries to be able to
% compute the correct near-boundary derivatives.
function dfieldovdn = diffxy(fieldmtx, fieldtype, BClist, ddir)
%   fieldmtx   - input field (2D matrix)
%   fieldtype  - type of input field: 'Ex','Ey','Ez','Hx','Hy' or 'Hz'
%   BClist     - boundary condition vector for the 2D field distribution
%                matrix, [left right bottom top]
%   ddir       - direction of derivative, n: 0=x, 1=y
%   dfieldovdx - output d/dx of field
PEC = 0; PMCh = 1; PBC = 2; PMC = 3; PECh = 4;          % [MP-21nov2003] Define boundary condition identifiers
BCl = BClist(1); BCr = BClist(2); BCd = BClist(3); BCu = BClist(4);

switch ddir
    case 0      % x-derivative
        switch lower(fieldtype)
            case {'ex','hy','hz'}
                if(BCl == PMC),  fieldmtx = [-fieldmtx(1,:); fieldmtx];                         % Add left field point for gradient only for PMC, PMCh or PBC
                    elseif(BCl == PMCh), fieldmtx = [zeros(1,size(fieldmtx,2)); fieldmtx];
                    elseif(BCl == PBC),  fieldmtx = fieldmtx([end 1:end],:);
                end
                if(BCr == PMC),  fieldmtx = [fieldmtx; -fieldmtx(end,:)];                       % Add right field point for gradient only for PMC and PMCh
                    elseif(BCr == PMCh), fieldmtx = [fieldmtx; zeros(1,size(fieldmtx,2))];
                end

            case {'ey','hx','ez'}
                if(BCl == PECh),  fieldmtx = [-fieldmtx(1,:); fieldmtx];                        % Add left field point for gradient only for PMC, PMCh or PBC
                    elseif(BCl == PEC),  fieldmtx = [zeros(1,size(fieldmtx,2)); fieldmtx];
                end
                if(BCr == PECh),  fieldmtx = [fieldmtx; -fieldmtx(end,:)];                      % Add right field point for gradient only for PMC and PMCh
                    elseif(BCr == PEC),  fieldmtx = [fieldmtx; zeros(1,size(fieldmtx,2))];
                    elseif(BCr == PBC),  fieldmtx = fieldmtx([1:end 1],:);
                end


            otherwise
                error([mfilename ' - diffxy(): unrecognized field type ' fieldtype]);
        end

    case 1      % y-derivative
        switch lower(fieldtype)
            case {'ex','hy','ez'}
                if(BCd == PECh),  fieldmtx = [-fieldmtx(:,1), fieldmtx];                        % Add bottom field point for gradient only for PMC, PMCh or PBC
                    elseif(BCd == PEC),  fieldmtx = [zeros(size(fieldmtx,1),1), fieldmtx];
                end
                if(BCu == PECh),  fieldmtx = [fieldmtx, -fieldmtx(:,end)];                      % Add top field point for gradient only for PMC and PMCh
                    elseif(BCu == PEC),  fieldmtx = [fieldmtx, zeros(size(fieldmtx,1),1)];
                    elseif(BCu == PBC),  fieldmtx = fieldmtx(:,[1:end 1]);
                end


            case {'ey','hx','hz'}
                if(BCd == PMC),  fieldmtx = [-fieldmtx(:,1), fieldmtx];                         % Add bottom field point for gradient only for PMC, PMCh or PBC
                    elseif(BCd == PMCh), fieldmtx = [zeros(size(fieldmtx,1),1), fieldmtx];
                    elseif(BCd == PBC),  fieldmtx = fieldmtx(:,[end 1:end]);
                end
                if(BCu == PMC),  fieldmtx = [fieldmtx, -fieldmtx(:,end)];                       % Add top field point for gradient only for PMC and PMCh
                    elseif(BCu == PMCh), fieldmtx = [fieldmtx, zeros(size(fieldmtx,1),1)];
                end

            otherwise
                error([mfilename ' - diffxy(): unrecognized field type ' fieldtype]);
        end

    otherwise
        error([mfilename ' - diffxy(): unrecognized derivative direction ' num2str(ddir)]);
end

dfieldovdn = diff(fieldmtx,1,1+ddir);  % d/dx fieldmtx or d/dy fieldmtx
%---------.---------.---------.---------.---------.---------.---------.--------|
