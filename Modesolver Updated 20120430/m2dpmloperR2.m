%---------.---------.---------.---------.---------.---------.---------.--------|
% Code:      2D-PML Cylindrical/Cartesian Vector-Field Modesolver for
%            Beta or Omega Eigenvalue
% Component: MATRIX OPERATOR GENERATOR (CORE ENGINE)
% Date:      June 4, 2002
% Author:    Milos Popovic
%
% Syntax:    [H] = m2dpmloperR2(er, Pr, Pz, rho, dr, dz, keig, OPTS)
%
%            Field grid is M x N (Er is (M-1)xN, Ez is Mx(N-1)).
%            Dielectric er is (2M-1)x(2N-1) as are Pr, Pz and rho.
%            dr, dz are scalar.
%
% Legend:
%
% OPTS.fieldmode = {'V'},'M','S' for vector, semi- and scalar modes
% OPTS.eigmode   = 'w',{'b'}     for omega or gamma (beta-r) eigenvalue mode
% OPTS.BC        = [left right down up] boundary conditions, where each can be:
%                   PEC PMCh PBC or PMC (regular = full step, h = halfstep, and
%                   acronyms refer to Perfect Electric and Magnetic and
%                   Periodic Boundary Conditions).
%                   note: PEC = 0, PMCh = 1, PBC = 2, PMC = 3.
% OPTS.TWOD      = Force 2D mode? {0} = 3D normal mode, 1 = 2D mode (TM)

% To do:
% ------
% - need to test and verify isolated Hrz and Hzr single operator pixels
%   that occur only when both BCleft/right == PBC *and* BCbot/top == PBC.
%   The current boundary condition test code passes regardless of whether
%   these pixels are turned on or zero; thus a targeted test is needed.
%
% Code Updates:
% -------------
% Jun  3, 2002 - [m2dpmloper.m] First MATLAB implementtn of vector-radial code.
%              - Added keig value that gets represents the *known* (w or beta)
%              - Added OPTS structure for passing in settings:
%                  OPTS.fmode = {'V'},'M','S' for vector, semi- and scalar mode.
%                  OPTS.emode = 'w',{'b'} for omega^2 or (beta*r)^2 eigenvalue.
% Jun 27, 2002 - Correction in operator Hrz and Hzz matrix elements; the CURL
%                operator was incorrectly derived for cylindrical coordinates.
% Jun 28, 2002 - Generalized dielectric profile to accept isotropic or diagonal
%                tensor (currently only diagonal tensor is supported).
%                er(:,:,1) is isotropic; er(:,:,1:3) gives xx,yy,zz components.
% Feb 13, 2003 - Multiple corrections to the Hrr, Hzz and a few to
%                off-diagonal elements. CORRECTED: boundary condition on
%                left/right for Hzz (spotted by Ben Williams); also
%                corrected normal fields to have a Neumann, not Dirichlet,
%                boundary condition, so that true PEC walls are used. This
%                is important for using the walls for symmetry but shouldn't
%                affect previous calculations with mode in middle and walls
%                far away.
% Feb 14, 2003 - Corrected previously uncorrected error in Hrz: the rho,
%                1/rho was to be removed from one of the terms in the Hrz
%                operator! (error was detected on June 27, 2002, but this
%                term was not corrected properly then).
%              - Adding OPTS.BC to specify which boundary condition is to
%                be used on each boundary of the computational domain.
%              - First working BC is PEC; working next on PMCh
%              - Renamed OPTS.fmode and OPTS.emode to fieldmode and
%                eigmode, to be compatible with m2wcyl and other code.
% Nov 20, 2003 - Work appears to have been in progress on PMCh BC on Feb
%                14-19, 2003. We'll continue to work on this now.
%              - Removed strange zeros added on Feb 18 - don't seem to be
%                needed! (they might be if the rows/columns corresponding
%                to zero fields weren't removed from the H matrix?)
% Nov 21, 2003 - PMCh BC seems to work for Hrr, Hzz, field matrix indexing
% Nov 22, 2003 - OPTS.TWOD option added; 0 = default (3D), 1 = 2D (erases
%                some matrix element components to force delT dot DT = 0).
% Nov 23, 2003 - MAJOR FIX: in Hrz, rho was supposed to be removed in two
%                places in the first term on June 27, 2002 when the
%                correction was found, but somehow the change was never
%                made. Thus, it was corrected NOW. Second, minor change: in
%                the second term, the rho_rz was absorbed into the
%                temporary variable C to make it more readable and save
%                flops (fewer multiplies).
%              - The code appears to maybe NOT NEED any corrections for
%                PMCh BC in Hrz and Hzr elements IF the zero element
%                rows/columns are erased from the matrix! Seems to
%                automatically work (makes sense looking at Yee grid).
% ----R2-------------
% Nov  2, 2006 - STARTING COMPLETE REWRITE to have full-step and half-step
%                PEC and PMC BCs (and potentially PBC); and normalize
%                dimensions to wavelength.  Tangential fields are added on
%                the edges of the comp domain on all sides.
%              - Added PECh BC and its macro definition (PECh = 4)
%              - Changing PMCh default boundary to be 1/2 pixel *outside*
%                domain, not inside domain as before. This is done to make
%                minimum variation of non-zero field points with varying
%                BCs - only the outsidemost tangential fields toggle - they
%                are on for PMC, (the new) PMCh, and off for PEC and PECh.
%              - PBCs are defined: the outermost tangential fields are on
%                on the left and bottom, and off on the right and top.
% Nov  3, 2006 - FINISHED new Hrr, Hzz and inhomogeneous term, and operator
%                construction and removal of unused field points; still
%                need Hrz and Hzr.
% Nov  4, 2006 - PBC works generally (tested)
%              - ***NOTE*** PBC does not work if index discontinuity is
%                right at domain boundary (this is likely because the
%                leftmost and rightmost pixel must be the same for PBC but
%                aren't in input - not sure about this though, since
%                rightmost pixel dielectric value should never be used with
%                PBC).
% Nov  6, 2006 - PEC, PECh, PMC, PMCh and PBC verified to work correctly
%                in general (except caveat for PBC given in Nov 4, '06 note
%                above) using boundary condition test script
%                's_run_test_boundary_conditions.m'.  It still needs to be
%                verified, because at very coarse discretizations the
%                relative error of some BCs rises from 1e-16 to 1e-11..
%                which could be just numerical error or could be a small
%                bug somewhere.
% Nov  8, 2006 - Everything works - off-diagonal operator blocks Hrz, Hzr
%                checked again for correctness.  All boundary condition
%                tests are passed.
% Apr 30, 2012 - Changed handling of bug-fixed spdiags function.  Inserted
%                a function spdiags2() into this code which is a wrapper
%                for spdiags() to correct its behavior.
%              - Removed "***NOTE***: requires a modified spdiags.m
%                function" in help text for this code, above.
%---------.---------.---------.---------.---------.---------.---------.--------|

function [H] = m2dpmloperR2(er, Pr, Pz, rho, dr, dz, keig, OPTS)
PEC = 0; PMCh = 1; PBC = 2; PMC = 3; PECh = 4;   % Define boundary condition identifiers
global VMODE; if isempty(VMODE), VMODE = 0; end  % Check diagnostic flag VMODE
if(VMODE>0), disp([mfilename ': R2 of vector-Helmholtz operator being used']); end
if nargin < 8                                    % Default mode: [MP] to fix w/ varargin
    OPTS.eigmode = 'b'; OPTS.fieldmode = 'V';    % beta, full-vector field mode
    OPTS.BC    = [PEC PEC PEC PEC];              % Default PEC boundary conditions
end
if ~isfield(OPTS,'BC'), OPTS.BC = [PEC PEC PEC PEC]; end  % Default BCs
if isempty(OPTS.BC),    OPTS.BC = [PEC PEC PEC PEC]; end
BCl = OPTS.BC(1); BCr = OPTS.BC(2); BCd = OPTS.BC(3); BCu = OPTS.BC(4); % Macros
if ~isfield(OPTS,'TWOD'), OPTS.TWOD = 0; end      % [MP-22Nov2003] Default TWOD = 0 (3D)
if isempty(OPTS.TWOD),    OPTS.TWOD = 0; end

t0 = clock;                                            % Start timing code.

% Resolve the "double-density" dielectric matrices into rr, ff and zz parts
M = (size(er,1)+1)/2; N = (size(er,2)+1)/2;            % Size of *basic* grid.
[er_rr, er_ff, er_zz]        = mresolve(er); clear er; % Split (2*M-1)x(2*N-1)
[Pr_rr, Pr_zr, Pr_zz, Pr_rz] = mresolve(Pr); clear Pr; % matrix into 3-4 small,
[Pz_rr, Pz_zr, Pz_zz, Pz_rz] = mresolve(Pz); clear Pz; % "interleaved" ones.
[rho_rr, rho_zr, rho_zz]     = mresolve(rho); clear rho;

% Pre-compute some handy values
dr2 = dr^2; dz2 = dz^2; drdz = dr*dz;

% Create operator components ("diagonals") in matrix form
% Hrr submatrix diagonals
%%TODO> - If boundary is PBC, make sure the edge dielectric constant is same
%%        on top and bottom (and/or left and right) / or print warning.
%%      - If left BC=PBC then right must be PBC also; same with bottom/top.
%%      - To test that the [safety] zero elements don't actually matter,
%%        set them to NaN, and then check the resulting operator matrix that
%%        there are no NaNs.  After check, these lines can be commented out.
%% NOTE: Lines below marked [safety] are supposed to be allowed to be removed (needs to be tested, but they were put there just to be safe)
%        Lines below marked [required] cannot be removed.
% Hrr element 1 - drdr part                                                % [MP] 2006-11-01 - recoded to have full-pixel and half-pixel PMC/PEC and PBC BCs
A = -1/dr2 * Pr_rr ./ (rho_rr.^2) * (~OPTS.TWOD);                          % [MP]Jun5-error fixed-rho "squared"; [Nov222003, OPTS.TWOD added]
B = rho_zr .* Pr_zr ./ er_ff;
C = rho_rr .* er_rr;
Hrr_l  =  A(:,1:N) .* B(1:M-1,1:N) .* C([M-1 1:M-2],1:N);                  % [MP] 2006-11-01 - recoded
%                                         ^-used only for PBC
Hrr_cl = -A(:,1:N) .* B(1:M-1,1:N) .* C(:,1:N);
Hrr_r  =  A(:,1:N) .* B(2:M,  1:N) .* C([2:M-1 1],1:N);
%                                              ^-used only for PBC
Hrr_cr = -A(:,1:N) .* B(2:M,  1:N) .* C(:,1:N);
clear A B C;
Hrr_lp = zeros(size(Hrr_l)); Hrr_rp = zeros(size(Hrr_r));                  % Extra left and right wrap-around matrix columns for PBCs only
switch(BCl)                                                                % [2006-11-02] COMPLETED
    case PBC
%        Hrr_lp(1,:) = Hrr_l(1,:); Hrr_l(1,:) = 0;                          % [required] Hrr_lp = new matrix diagonal for left PBC wrapping
        Hrr_lp(1,1:N-1) = Hrr_l(1,1:N-1); Hrr_l(1,:) = 0;                  % [required; j=N left zero for safety (that can be ignored as in line above)] Hrr_lp = new matrix diagonal for left PBC wrapping
    case {PMC,PMCh}                                                        % new PMCh case goes here
        Hrr_l(1,:) = 0;                                                    % [required]
        if(BCl==PMC), Hrr_cl(1,:) = 2*Hrr_cl(1,:); end                     % [required] for PMC needs doubling, for PMCh stays the same
    case {PEC,PECh}                                                        % old PMCh case was here
        Hrr_l(1,:) = 0; Hrr_cl(1,:) = 0;                                   % [required] to avoid wrapping to other side of domain
        if(BCl==PECh), Hrr_r(1,:) = 2*Hrr_r(1,:); Hrr_cr(1,:) = 2*Hrr_cr(1,:); end  % [required] for PECh needs this doubling (symmetry); for PEC keep original value
    otherwise
        error(['Unknown left boundary condition type ' num2str(BCl)]);
end
switch(BCr)                                                                % [2006-11-02] COMPLETED
    case PBC
%        Hrr_rp(M-1,:) = Hrr_r(M-1,:); Hrr_r(M-1,:) = 0;                    % [required] Hrr_rp = new matrix diagonal for right PBC wrapping
        Hrr_rp(M-1,1:N-1) = Hrr_r(M-1,1:N-1); Hrr_r(M-1,:) = 0;            % [required; j=N left zero for safety (that can be ignored as in line above)] Hrr_rp = new matrix diagonal for right PBC wrapping
    case {PMC,PMCh}                                                        % new PMCh case goes here
        Hrr_r(M-1,:) = 0;                                                  % [required]
        if(BCr==PMC), Hrr_cr(M-1,:) = 2*Hrr_cr(M-1,:); end                 % [required] for PMC needs doubling, for PMCh stays the same
    case {PEC,PECh}                                                        % old PMCh case was here
        Hrr_r(M-1,:) = 0; Hrr_cr(M-1,:) = 0;                               % [required] to avoid wrapping to other side of domain
        if(BCr==PECh), Hrr_l(M-1,:) = 2*Hrr_l(M-1,:); Hrr_cl(M-1,:) = 2*Hrr_cl(M-1,:); end  % [required] for PECh needs this doubling (symmetry); for PEC keep original value
    otherwise
        error(['Unknown right boundary condition type ' num2str(BCr)]);
end

% Hrr element 2 - dzdz part                                                % [MP] 2006-11-01 - recoded to have full-pixel and half-pixel PMC/PEC and PBC BCs
Hrr_u = -1/dz2 * Pz_rr(:,1:N) .* Pz_rz(:,[1:N-1 1]);  Hrr_cu = -Hrr_u;     % (M-1)x(N-1+1)
%                                               ^-filler (not used)
Hrr_d = -1/dz2 * Pz_rr(:,1:N) .* Pz_rz(:,[N-1 1:N-1]);  Hrr_cd = -Hrr_d;   % (M-1)x(N-1+1)
%                                          ^-not used in Hrr_d (filler); is used for PBC in Hrr_cd
switch(BCu)                                                                % [2006-11-02] COMPLETED
    case PBC
        Hrr_u(:,N) = 0; Hrr_cu(:,N) = 0; Hrr_cd(:,N) = 0; Hrr_d(:,N)  = 0; % [rows zeroed for safety] Edge Er field rows will be removed from operator, so this line isn't strictly necessary
        Hrr_l(:,N) = 0; Hrr_cl(:,N) = 0; Hrr_cr(:,N) = 0; Hrr_r(:,N)  = 0; % [safety]
        % [safety] - Could also set Hrr_u(:,N-1) = 0 for the main up
        % diagonal, but keep a new Hrr_up(:,N-1) (see left/right case) for
        % the added "wrapped" diagonal
    case {PMC,PMCh}                                                        % new PMCh case goes here
        Hrr_u(:,N) = 0;                                                    % [row N zeroed for safety] Reaches above last row, so will run off the matrix
        Hrr_cu(:,N) = 0;                                                   % [required]
        if(BCu==PMC), Hrr_d(:,N) = 2*Hrr_d(:,N); Hrr_cd(:,N) = 2*Hrr_cd(:,N); end  % [required] For Neumann full-pixel PMC BC double value (symmetry);  for PMCh keep unchanged.
    case {PEC,PECh}                                                        % old PMCh case was here
        Hrr_u(:,N) = 0; Hrr_cu(:,N) = 0; Hrr_cd(:,N) = 0; Hrr_d(:,N)  = 0; % [safety PEC & PECh] rows zeroed for safety; Edge Er field rows will be removed from operator, so this line isn't strictly necessary
        Hrr_l(:,N) = 0; Hrr_cl(:,N) = 0; Hrr_cr(:,N) = 0; Hrr_r(:,N)  = 0; % [safety PEC & PECh]
        Hrr_u(:,N-1) = 0;                                                  % [safety for PEC and PECh] these field points will be removed from operator so this shouldn't show up
        if(BCu==PECh), Hrr_cu(:,N-1) = 2*Hrr_cu(:,N-1); end                % [required] for PECh needs this doubling (symmetry); for PEC keep original value
    otherwise
        error(['Unknown upper boundary condition type ' num2str(BCu)]);
end
switch(BCd)                                                                % [2006-11-02] COMPLETED
    case PBC
        % Do nothing
    case {PMC,PMCh}                                                        % new PMCh case goes here
        Hrr_d(:,1) = 0;                                                    % [row 1 zero for safety] Bottom Er field row will be removed from operator, so zeroing first line isn't strictly necessary
        Hrr_cd(:,1) = 0;                                                   % [required]
        if(BCd==PMC), Hrr_u(:,1) = 2*Hrr_u(:,1); Hrr_cu(:,1) = 2*Hrr_cu(:,1); end  % [required] For Neumann full-pixel PMC BC double value (symmetry);  for PMCh keep unchanged.
    case {PEC,PECh}                                                        % old PMCh case was here
        Hrr_u(:,1) = 0; Hrr_cu(:,1) = 0; Hrr_cd(:,1) = 0; Hrr_d(:,1)  = 0; % [safety] Edge Er field rows will be removed from operator, so this line isn't strictly necessary
        Hrr_l(:,1) = 0; Hrr_cl(:,1) = 0; Hrr_cr(:,1) = 0; Hrr_r(:,1)  = 0; % [safety PEC and PECh]
        Hrr_d(:,2) = 0;                                                    % [safety for PEC and PECh] Bottom Er field row will be removed from operator, so zeroing first line isn't strictly necessary
        if(BCd==PECh), Hrr_cd(:,2) = 2*Hrr_cd(:,2); end                    % [required] for PECh needs this doubling (symmetry); for PEC keep original value
    otherwise
        error(['Unknown lower boundary condition type ' num2str(BCd)]);
end

Hrr_c = Hrr_cl + Hrr_cr + Hrr_cu + Hrr_cd;                                 % Construct center-field element of 5-point stencil block, from its 4 subparts
clear Hrr_cl Hrr_cr Hrr_cu Hrr_cd;

% Construct complete Hrr operator block
%Hrr = spdiags([Hrr_d(:) Hrr_l(:) Hrr_c(:) Hrr_r(:) Hrr_u(:)], ...         % Standard form without support for PBCs
%              [(M-1) 1 0 -1 -(M-1)], (M-1)*N, (M-1)*N).';
% For PBCs only -v                 --v                                  --v                --v
%Hrr = spdiags([Hrr_u(:) Hrr_d(:) Hrr_rp(:) Hrr_l(:) Hrr_c(:) Hrr_r(:) Hrr_lp(:) Hrr_u(:) Hrr_d(:)], ...     % Added 4 diagonals for PBCs
%              [(-(M-1)+(M-1)*N) (M-1) (M-2) 1 0 -1 -(M-2) -(M-1) ((M-1)-(M-1)*N)], (M-1)*N, (M-1)*N).';
Hrr = spdiags2([Hrr_u(:)*(BCu==PBC) Hrr_d(:) Hrr_rp(:) Hrr_l(:) Hrr_c(:) Hrr_r(:) Hrr_lp(:) Hrr_u(:) Hrr_d(:)*(BCd==PBC)], ...     % Added 4 diagonals for PBCs
              [(-(M-1)+(M-1)*(N-1)) (M-1)    (M-2)       1        0        -1    -(M-2)     -(M-1) (+(M-1)-(M-1)*(N-1))], (M-1)*N, (M-1)*N).';
clear Hrr_l Hrr_r Hrr_u Hrr_d Hrr_lp Hrr_rp;
%figure; spy(Hrr);                                                         % [MP] For debugging.
%----


% Hzz submatrix diagonals
% Hzz element 1 - drdr part                                                % [MP] 2006-11-01 - recoded to have full-pixel and half-pixel PMC/PEC and PBC BCs
A1 = Pr_zz/dr2 ./ rho_zz;                                                  % [MP] Jun 27 correction
Hzz_l = -A1(1:M,:) .* Pr_rz([M-1 1:M-1],:) .* rho_rr([M-1 1:M-1],1:N-1);  Hzz_cl = -Hzz_l; % [MP] Should be rho_rz!, but rho_rz is not obtained from mresolve(rho..) above.
%                             ^-used only for PBC------^
Hzz_r = -A1(1:M,:) .* Pr_rz([1:M-1 1],:)   .* rho_rr([1:M-1 1],1:N-1);  Hzz_cr = -Hzz_r;
%                                  ^-used only for PBC------^
clear A;
% Hzz element 2 - dzdz part                                                % [MP] 2006-11-01 - recoded to have full-pixel and half-pixel PMC/PEC and PBC BCs
B = -Pz_zz/dz2  * (~OPTS.TWOD); C = Pz_zr ./ er_ff;                        % [MP-Nov222003] OPTS.TWOD added
Hzz_u  =  B(1:M,:) .* C(1:M,2:N)   .* er_zz(1:M,[2:N-1 1]);
%                                                      ^-used only for PBC
Hzz_cu = -B(1:M,:) .* C(1:M,2:N)   .* er_zz(1:M,:);
Hzz_d  =  B(1:M,:) .* C(1:M,1:N-1) .* er_zz(1:M,[N-1 1:N-2]);
%                                                 ^-used only for PBC
Hzz_cd = -B(1:M,:) .* C(1:M,1:N-1) .* er_zz(1:M,:);
clear B C;

Hzz_lp = zeros(size(Hzz_l)); Hzz_rp = zeros(size(Hzz_r));                  % Extra left and right wrap-around matrix columns for PBCs only
switch(BCl)                                                                % [2006-11-02] COMPLETED
    case PBC
        Hzz_lp(1,:) = Hzz_l(1,:); Hzz_l(1,:) = 0;                          % [required] Hzz_lp = new matrix diagonal for left PBC wrapping
    case {PMC,PMCh}
        Hzz_l(1,:) = 0;                                                    % [row N zeroed for safety] Top Er field row will be removed from operator, so this line isn't strictly necessary
        Hzz_cl(1,:) = 0;                                                   % [required]
        if(BCl==PMC), Hzz_r(1,:) = 2*Hzz_r(1,:); Hzz_cr(1,:) = 2*Hzz_cr(1,:); end  % [required]
    case {PEC,PECh}
        Hzz_l(1,:) = 0; Hzz_cl(1,:) = 0; Hzz_cr(1,:) = 0; Hzz_r(1,:)  = 0; % [safety]
        Hzz_u(1,:) = 0; Hzz_cu(1,:) = 0; Hzz_cd(1,:) = 0; Hzz_d(1,:)  = 0; % [safety PEC and PECh]
        Hzz_l(2,:) = 0;                                                    % [safety for PEC and PECh]
        if(BCl==PECh), Hzz_cl(2,:) = 2*Hzz_cl(2,:); end                    % [required] for PECh needs this doubling (symmetry); for PEC keep original value
    otherwise
        error(['Unknown left boundary condition type ' num2str(BCl)]);
end
switch(BCr)                                                                % [2006-11-02] COMPLETED
    case PBC
        Hzz_l(M,:) = 0; Hzz_cl(M,:) = 0; Hzz_cr(M,:) = 0; Hzz_r(M,:)  = 0; % [safety]
        Hzz_u(M,:) = 0; Hzz_cu(M,:) = 0; Hzz_cd(M,:) = 0; Hzz_d(M,:)  = 0; % [safety]
        Hzz_rp(M-1,:) = Hzz_r(M-1,:); Hzz_r(M-1,:) = 0;                    % [required] Hzz_rp = new matrix diagonal for right PBC wrapping
    case {PMC,PMCh}
        Hzz_r(M,:) = 0;                                                    % [row N zeroed for safety] Top Er field row will be removed from operator, so this line isn't strictly necessary
        Hzz_cr(M,:) = 0;                                                   % [required]
        if(BCr==PMC), Hzz_l(M,:) = 2*Hzz_l(M,:); Hzz_cl(M,:) = 2*Hzz_cl(M,:); end  % [required] For Neumann full-pixel PMC BC double value (symmetry);  for PMCh keep unchanged.
    case {PEC,PECh}
        Hzz_l(M,:) = 0; Hzz_cl(M,:) = 0; Hzz_cr(M,:) = 0; Hzz_r(M,:)  = 0; % [safety PEC & PECh]
        Hzz_u(M,:) = 0; Hzz_cu(M,:) = 0; Hzz_cd(M,:) = 0; Hzz_d(M,:)  = 0; % [safety PEC & PECh]
        Hzz_r(M-1,:) = 0;                                                  % [safety for PEC and PECh] these field points will be removed from operator so this shouldn't show up
        if(BCr==PECh), Hzz_cr(M-1,:) = 2*Hzz_cr(M-1,:); end                % [required] for PECh needs this doubling (symmetry); for PEC keep original value
    otherwise
        error(['Unknown right boundary condition type ' num2str(BCr)]);
end
switch(BCu)                                                                % [2006-11-02] COMPLETED
    case PBC
        % Do nothing
    case {PMC,PMCh}
        Hzz_u(:,N-1) = 0;                                                  % [required]
        if(BCu==PMC), Hzz_cu(:,N-1) = 2*Hzz_cu(:,N-1); end                 % [required] for PMC needs doubling, for PMCh stays the same
    case {PEC,PECh}
        Hzz_u(:,N-1) = 0; Hzz_cu(:,N-1) = 0;                               % [required] to avoid wrapping to other side of domain
        if(BCu==PECh), Hzz_d(:,N-1) = 2*Hzz_d(:,N-1); Hzz_cd(:,N-1) = 2*Hzz_cd(:,N-1); end  % [required] for PECh needs this doubling (symmetry); for PEC keep original value
    otherwise
        error(['Unknown upper boundary condition type ' num2str(BCu)]);
end
switch(BCd)                                                                % [2006-11-02] COMPLETED
    case PBC
        % Do nothing
    case {PMC,PMCh}
        Hzz_d(:,1) = 0;                                                    % [required]
        if(BCd==PMC), Hzz_cd(:,1) = 2*Hzz_cd(:,1); end                     % [required] for PMC needs doubling, for PMCh stays the same
    case {PEC,PECh}
        Hzz_d(:,1) = 0; Hzz_cd(:,1) = 0;                                   % [required] to avoid wrapping to other side of domain
        if(BCd==PECh), Hzz_u(:,1) = 2*Hzz_u(:,1); Hzz_cu(:,1) = 2*Hzz_cu(:,1); end  % [required] for PECh needs this doubling (symmetry); for PEC keep original value
    otherwise
        error(['Unknown lower boundary condition type ' num2str(BCd)]);
end

Hzz_c = Hzz_cl + Hzz_cr + Hzz_cu + Hzz_cd;                                 % Construct center-field element of 5-point stencil block, from its 4 subparts
clear Hzz_cl Hzz_cr Hzz_cu Hzz_cd;

% Construct complete Hzz operator block
%Hzz = spdiags([Hzz_d(:) Hzz_l(:) Hzz_c(:) Hzz_r(:) Hzz_u(:)], ...         % Standard form without support for PBCs
%              [M 1 0 -1 -M], M*(N-1), M*(N-1)).';
% NOTE: it is not clear if the Hzz_u and Hzz_d for PBC BCs need to be
% set to zero for non-PBC BCs - it appears ok without, but needs to be
% checked.
% Shown below is the (BCu==PBC) flag that was required in the Hrr case.
%Hrr = spdiags([Hrr_u(:)*(BCu==PBC) Hrr_d(:) Hrr_rp(:) Hrr_l(:) Hrr_c(:) Hrr_r(:) Hrr_lp(:) Hrr_u(:) Hrr_d(:)*(BCd==PBC)], ...     % Added 4 diagonals for PBCs
%              [(-(M-1)+(M-1)*(N-1)) (M-1)    (M-2)       1        0        -1    -(M-2)     -(M-1) (+(M-1)-(M-1)*(N-1))], (M-1)*N, (M-1)*N).';
% For PBCs only -v                    --v                                  --v                  --v
Hzz = spdiags2([Hzz_u(:)     Hzz_d(:) Hzz_rp(:) Hzz_l(:) Hzz_c(:) Hzz_r(:) Hzz_lp(:) Hzz_u(:)   Hzz_d(:)], ...     % Added 4 diagonals for PBCs
              [(-M+M*(N-1))    M      (M-2)       1        0        -1      -(M-2)    -M    (M-M*(N-1))], M*(N-1), M*(N-1)).';
clear Hzz_l Hzz_r Hzz_u Hzz_d Hzz_lp Hzz_rp;
%figure; spy(Hzz);                                                         % [MP] For debugging.
%----



% Off-diagonal submatrices only for full-vector solutions:
if(OPTS.fieldmode == 'V') % 'V' = Vector, 'M' = Semi-Vector (quasi-...), 'S' = Scalar
    % Hrz submatrix diagonals
%% [OLD CODE - Deleted another commented out version of below operator forming code - see previous m2dpmloper.m versions]
%% NEW CODE (corrected and simplified) - 23 nov 2003 %%
%     A1 = Pz_rr(:,2:N-1)/drdz; A2 = A1 .* Pr_rz(:,2:N-1);                    % [MP-23nov2003] MAJOR FIX! Removed rho's in first term of Hrz. This was claimed done on June 27, 2002 but wasn't!!
%     B  = Pr_rr(:,2:N-1)/drdz ./ (rho_rr(:,2:N-1).^2);                       % Can use variable A from part for Hrr
%     C  = rho_zr(:,2:N-1).^2 .* Pz_zr(:,2:N-1) ./ er_ff(:,2:N-1);            % [MP-23nov2003] absorbed rho_zr into C to simplify lines below
%     Hrz_c = -[zeros(1,N-2); A2(2:M-1,:)] + ...                              % [MP-23nov2003] MAJOR FIX!! Rho removed from Hrz 1st term
%              B  .* [zeros(1,N-2); C(2:M-1,:)] .* er_zz(1:M-1,2:N-1) * (~OPTS.TWOD); % [MP22nov2003-OPTS.2D], [23nov03] rho_rz absorbed into C
%     Hrz_r = +[A2(1:M-2,:); zeros(1,N-2)] - ...                              % [MP-23nov2003] MAJOR FIX!! Rho removed from Hrz 1st term
%              B  .* [C(2:M-1,:); zeros(1,N-2)] .* er_zz(2:M,2:N-1) * (~OPTS.TWOD);   % [MP22nov2003-OPTS.2D], [23nov03] rho_rz absorbed into C
%     A2 = A1 .* Pr_rz(:,1:N-2);
%     clear A1;
%     Hrz_d = +[zeros(1,N-2); A2(2:M-1,:)] - ...                              % or rho_zr(..,2:N-1)    % [MP-23nov2003] MAJOR FIX!! Rho removed from Hrz 1st term
%              B  .* [zeros(1,N-2); C(2:M-1,:)] .* er_zz(1:M-1,1:N-2) * (~OPTS.TWOD); % [MP22nov2003-OPTS.2D], [23nov03] rho_rz absorbed into C
%     Hrz_rd= -[A2(1:M-2,:); zeros(1,N-2)] + ...                              % [MP-23nov2003] MAJOR FIX!! Rho removed from Hrz 1st term
%              B  .* [C(2:M-1,:); zeros(1,N-2)] .* er_zz(2:M,1:N-2) * (~OPTS.TWOD);   % [MP22nov2003-OPTS.2D], [23nov03] rho_rz absorbed into C
%     clear A2 B C;
%% [NEW CODE PMCh.. - Deleted another commented out version of above operator forming code - see previous m2dpmloper.m versions]

% [2006-11-07] New code for Hrz operator block
    A1 = Pz_rr(:,1:N)/drdz;
    A2 = A1 .* Pr_rz(:,[1:N-1 1]);
%                             ^-dummy(unused)
    B  = Pr_rr(:,1:N)/drdz ./ (rho_rr(:,1:N).^2);                           % TODO: Can use variable A from part for Hrr?
    C  = rho_zr(:,1:N).^2 .* Pz_zr(:,1:N) ./ er_ff(:,1:N);
    Hrz_lu = - A2 + B .* C(1:M-1,1:N) .* er_zz(1:M-1,[1:N-1 1]);            % Left-up stencil quarter
    Hrz_ru = + A2 - B .* C(2:M,  1:N) .* er_zz(2:M,  [1:N-1 1]);
%                                                           ^-dummy(unused)
    A2 = A1 .* Pr_rz(:,[N-1 1:N-1]);
%                        ^-used for PBC only
    clear A1;
    Hrz_ld = + A2 - B .* C(1:M-1,1:N) .* er_zz(1:M-1,[N-1 1:N-1]);
    Hrz_rd = - A2 + B .* C(2:M,  1:N) .* er_zz(2:M,  [N-1 1:N-1]);
%                                                      ^-used for PBC only
    clear A2 B C;
% * (~OPTS.TWOD);

    %--- Begin set Hrz boundary conditions
    Hrz_rdp = zeros(size(Hrz_rd)); Hrz_rup = zeros(size(Hrz_ru));           % Extra right wrap-around matrix columns for PBCs only
    switch(BCl)                                                             % [2006-11-07] COMPLETED
        case {PBC,PMC,PMCh}
            % Do nothing
        case {PEC,PECh}
            Hrz_lu(1,:) = 0; Hrz_ld(1,:) = 0;                               % [safety]
            if(BCl==PECh), Hrz_ru(1,:) = 2*Hrz_ru(1,:); Hrz_rd(1,:) = 2*Hrz_rd(1,:); end % [required] for PECh needs this doubling (symmetry); for PEC keep original value
        otherwise
            error(['Unknown left boundary condition type ' num2str(BCl)]);
    end
    switch(BCr)                                                             % [2006-11-07] COMPLETED
        case PBC
            Hrz_rdp(M-1,:) = Hrz_rd(M-1,:); Hrz_rd(M-1,:) = 0;              % [required] Hrz_rdp and Hrz_rup = new matrix diagonals for right PBC wrapping
            Hrz_rup(M-1,:) = Hrz_ru(M-1,:); Hrz_ru(M-1,:) = 0;              % [required]
        case {PMC,PMCh}
            % Do nothing
        case {PEC,PECh}
            Hrz_ru(M-1,:) = 0; Hrz_rd(M-1,:) = 0;                           % [safety]
            if(BCr==PECh), Hrz_lu(M-1,:) = 2*Hrz_lu(M-1,:); Hrz_ld(M-1,:) = 2*Hrz_ld(M-1,:); end % [required] for PECh needs this doubling (symmetry); for PEC keep original value
        otherwise
            error(['Unknown right boundary condition type ' num2str(BCr)]);
    end
    switch(BCu)                                                             % [2006-11-07] COMPLETED
        case {PMC,PMCh}
            Hrz_lu(:,N) = 0; Hrz_ru(:,N) = 0;                               % [safety]
            if(BCu==PMC), Hrz_ld(:,N) = 2*Hrz_ld(:,N); Hrz_rd(:,N) = 2*Hrz_rd(:,N); end % [required] for PMC needs doubling, for PMCh stays the same
        case {PBC,PEC,PECh}
            Hrz_lu(:,N) = 0; Hrz_ru(:,N) = 0; Hrz_ld(:,N) = 0; Hrz_rd(:,N) = 0; % [safety] top row of tangential E-fields will be auto removed from matrix
            % Nothing special to do for PECh!
            % Do nothing for PBC either.
        otherwise
            error(['Unknown upper boundary condition type ' num2str(BCu)]);
    end
    switch(BCd)                                                             % [2006-11-07] COMPLETED
        case PBC
            % Do nothing - weighting factors for fields below boundary must come from top of domain
        case {PMC,PMCh}
            Hrz_ld(:,1) = 0; Hrz_rd(:,1) = 0;                               % [safety] bottom row of tangential E-fields will be auto removed from matrix
            if(BCd==PMC), Hrz_lu(:,1) = 2*Hrz_lu(:,1); Hrz_ru(:,1) = 2*Hrz_ru(:,1); end % [required] for PMC needs doubling, for PMCh stays the same
        case {PEC,PECh}
            Hrz_lu(:,1) = 0; Hrz_ru(:,1) = 0; Hrz_ld(:,1) = 0; Hrz_rd(:,1) = 0; % [safety] bottom row of tangential E-fields will be auto removed from matrix
        otherwise
            error(['Unknown lower boundary condition type ' num2str(BCd)]);
    end
    %--- End of setting Hrz boundary conditions

% Not needed:
    Hrz_lu = [Hrz_lu; zeros(1,N)]; Hrz_ru = [Hrz_ru; zeros(1,N)]; Hrz_rup = [Hrz_rup; zeros(1,N)];          % Fill left/right with zeros to match the input Ez field size before diagonalizing block
    Hrz_ld = [Hrz_ld; zeros(1,N)]; Hrz_rd = [Hrz_rd; zeros(1,N)]; Hrz_rdp = [Hrz_rdp; zeros(1,N)];


%%    Hrz = spdiags([Hrz_d(:) Hrz_rd(:) Hrz_c(:) Hrz_r(:)], ...
%%                  [1 0 1-(M-2+1) 0-(M-2+1)], (M-2)*(N-1) + (N-2), (M-1)*(N-2)).';
%%    clear Hrz_c Hrz_r Hrz_d Hrz_rd;
%%    Hrz = spdiags([Hrz_ld(:) Hrz_rd(:) Hrz_lu(:) Hrz_ru(:)], ...
%%                  [   1         0      1-(M-2+1) 0-(M-2+1)], (M-2)*(N-1) + (N-2), (M-1)*(N-2)).';
%%%                                                            input size -^             ^- output size
%%    Hrz = spdiags([Hrz_ld(:) Hrz_rd(:) Hrz_lu(:) Hrz_ru(:)], ...
%%                  [   0         -1        -M      -(M+1)    ], M*(N+1), (M-1)*N).';
%%%                                                            input size -^             ^- output size
%    Hrz = spdiags([Hrz_rdp(:)*(BCr==PBC) Hrz_ld(:) Hrz_rd(:) Hrz_rup(:)*(BCr==PBC) Hrz_lu(:) Hrz_ru(:) Hrz_ld(:)*(BCd==PBC)   Hrz_rd(:)*(BCd==PBC) ], ...  % Complete BC support in Hrz
%                  [   M-2                   0         -1        -2                    -M      -(M+1)   -(N-2)*M             -((N-2)*M+1)           ], M*(N+1), (M-1)*N).';
%%                                                            input size -^             ^- output size
    Hrz = spdiags2([Hrz_rdp(:)*(BCr==PBC) Hrz_ld(:) Hrz_rd(:) Hrz_rup(:)*(BCr==PBC) Hrz_lu(:) Hrz_ru(:) Hrz_ld(:)*(BCd==PBC)   Hrz_rd(:)*(BCd==PBC) ], ...  % Complete BC support in Hrz
                  [  (M+M-2)                 M        M-1       M-2                    0       -1          -(N-2)*M             -((N-2)*M+1)       ], ...
                  M*(N-1), M*N).'; % M*N to be cut to (M-1)*N by row erasing
%        input size -^      ^- output size is (M-1)*N, needs reduction by row erasing
    Hrz(M-1,(N-2)*M+1) = Hrz_rd(M-1,1)*(BCr==PBC)*(BCd==PBC); % Single matrix element for bottom right corner wrapping to top left corner when both PBCs are on.
    clear Hrz_ld Hrz_lu Hrz_rd Hrz_ru Hrz_rdp Hrz_rup;
%    nidx = [];                                                              % Filtering out blank lines in Hrz to make the submatrix properly
%%     for m = 1:N-1                                                           % "non-diagonal" (since Er is larger in r than Ez by 1 pixel).
%%         nidx = [nidx (m-1)*((M-2)+1)+(1:M-2)];
%%     end
%    nidx = M+[1:M*(N-1)];   % Remove top and bottom row immediately (this is done here because the bottom row was needed to get the second PBC diagonal)
%    Hrz = Hrz(:,nidx);      % [MP] Reformat the Hrz matrix eliminating the above empty lines.
    nidx = setdiff((1:M*N), M*(1:N));
    Hrz = Hrz(nidx,:);        % [MP] Reformat the Hrz matrix eliminating the above empty lines.
%    spy(Hrz);

    % Hzr submatrix diagonals
%     A1 = Pr_zz(2:M-1,:)/drdz ./ rho_zz(2:M-1,:); A2 = A1 .* Pz_rz(2:M-1,:);
%     B1 = Pz_zz(2:M-1,:)/drdz ./ rho_zz(2:M-1,:);
%     B2 = B1 .* Pr_zr(2:M-1,1:N-1) ./ er_ff(2:M-1,1:N-1);
%     Hzr_c = -A2(:,2:N-1) .* rho_rr(2:M-1,2:N-1) ...
%             +B2(:,2:N-1) .* rho_rr(2:M-1,2:N-1) .* er_rr(2:M-1,2:N-1) * (~OPTS.TWOD); % [MP22nov2003]
%     A1 = A1 .* Pz_rz(1:M-2,:);
%     Hzr_l =  A1(:,2:N-1) .* rho_rr(1:M-2,2:N-1) ...
%             -B2(:,2:N-1) .* rho_rr(1:M-2,2:N-1) .* er_rr(1:M-2,2:N-1) * (~OPTS.TWOD); % [MP22nov2003]
%     clear B2; B1 = B1 .* Pr_zr(2:M-1,2:N) ./ er_ff(2:M-1,2:N);
%     Hzr_u =  A2(:,1:N-2) .* rho_rr(2:M-1,2:N-1) ...
%             -B1(:,1:N-2) .* rho_rr(2:M-1,2:N-1) .* er_rr(2:M-1,2:N-1) * (~OPTS.TWOD); % [MP22nov2003]
%     clear A2;
%     Hzr_lu= -A1(:,1:N-2) .* rho_rr(1:M-2,2:N-1) ...
%             +B1(:,1:N-2) .* rho_rr(1:M-2,2:N-1) .* er_rr(1:M-2,2:N-1) * (~OPTS.TWOD); % [MP22nov2003]
%     clear A1 B1;
% 
%     A1 = zeros(1,size(Hzr_l,2));
%     Hzr_l = [Hzr_l; A1]; Hzr_c = [Hzr_c; A1]; Hzr_lu = [Hzr_lu; A1]; Hzr_u = [Hzr_u; A1];
%     Hzr = spdiags([[Hzr_l(:);0] [0;Hzr_c(:)] [Hzr_lu(:);0] [0;Hzr_u(:)]], ...
%                   [-(M-2)-1 -(M-2)+1-1 0 1], (M-2)*(N-1)+(N-2), (M-1)*(N-2));
%     clear Hzr_c Hzr_l Hzr_u Hzr_lu;
%     nidx = [];      % Filtering out blank horiz lines in Hzr to get proper submatrix
%     for m = 1:N-1   % "non-diagonal" (since Er is larger in r than Ez by 1 pixel).
%         nidx = [nidx (m-1)*((M-2)+1)+(1:M-2)];
%     end
%     Hzr = Hzr(nidx,:);

% [2006-11-07] New code for Hzr operator block
%%%if(0)
%%%    Hzr = Hrz';
%%%else
    A1 = Pr_zz(1:M,:)/drdz ./ rho_zz(1:M,:); A2 = A1 .* Pz_rz([1:M-1 1],:) .* rho_rr([1:M-1 1],1:N-1);
%                                                                    ^----dummy (unused)----^
    B1 = Pz_zz(1:M,:)/drdz ./ rho_zz(1:M,:);                                % z-dimension is (N-1) in length
    B2 = B1 .* Pr_zr(1:M,1:N-1) ./ er_ff(1:M,1:N-1);
    C  = rho_rr .* er_rr;

    Hzr_rd = -A2 + B2 .* C([1:M-1 1],1:N-1); % * (~OPTS.TWOD);              % Substituting rho_rr for the non-existent rho_rz (since it should have no z-dependence)
%                                 ^-dummy (unused)
    A1 = A1 .* Pz_rz([M-1 1:M-1],:) .* rho_rr([M-1 1:M-1],1:N-1);
%                      ^-used for PBC only -----^
    Hzr_ld = +A1 - B2 .* C([M-1 1:M-1],1:N-1); % * (~OPTS.TWOD);
%                            ^-used for PBC only
    clear B2; B1 = B1 .* Pr_zr(1:M,2:N) ./ er_ff(1:M,2:N);
    Hzr_ru = +A2 - B1 .* C([1:M-1 1],2:N); % * (~OPTS.TWOD);
%                                 ^-dummy (unused)
    clear A2;
    Hzr_lu = -A1 + B1 .* C([M-1 1:M-1],2:N); % * (~OPTS.TWOD);
%                            ^-used for PBC only
    clear A1 B1 C;

    %--- Begin set Hzr boundary conditions
    Hzr_ldp = zeros(size(Hzr_ld)); Hzr_lup = zeros(size(Hzr_lu));           % Extra left wrap-around matrix columns for PBCs only
    switch(BCd)                                                             % [2006-11-07] COMPLETED
        case {PBC,PMC,PMCh}
            % Do nothing
        case {PEC,PECh}
            Hzr_ld(:,1) = 0; Hzr_rd(:,1) = 0;                               % [safety] should be wiped automatically since the input field they lie on is zero
            if(BCd==PECh), Hzr_lu(:,1) = 2*Hzr_lu(:,1); Hzr_ru(:,1) = 2*Hzr_ru(:,1); end % [required] for PECh needs this doubling (symmetry); for PEC keep original value
        otherwise
            error(['Unknown lower boundary condition type ' num2str(BCd)]);
    end
    switch(BCu)                                                             % [2006-11-07] COMPLETED
        case PBC
            % Do nothing - top tangential E-field factors run off the
            % matrix "up-diagonal", but will show up in the secondary
            % diagonal "up-prime" used only for the PBC
        case {PMC,PMCh}
            % Do nothing
        case {PEC,PECh}
            Hzr_lu(:,N-1) = 0; Hzr_ru(:,N-1) = 0;                           % [safety] since top row tangential E fields multiply a zero input field (which is removed later from the matrix)
            if(BCu==PECh), Hzr_ld(:,N-1) = 2*Hzr_ld(:,N-1); Hzr_rd(:,N-1) = 2*Hzr_rd(:,N-1); end % [required] for PECh needs this doubling (symmetry); for PEC keep original value
        otherwise
            error(['Unknown upper boundary condition type ' num2str(BCu)]);
    end
    switch(BCr)                                                             % [2006-11-07] COMPLETED
        case {PMC,PMCh}
            Hzr_ru(M,:) = 0; Hzr_rd(M,:) = 0;                               % [required]
            if(BCr==PMC), Hzr_lu(M,:) = 2*Hzr_lu(M,:); Hzr_ld(M,:) = 2*Hzr_ld(M,:); end % [required] for PMC needs doubling, for PMCh stays the same
        case {PBC,PEC,PECh}
            Hzr_lu(M,:) = 0; Hzr_ld(M,:) = 0; Hzr_ru(M,:) = 0; Hzr_rd(M,:) = 0; % [safety] top row of tangential E-fields will be auto removed from matrix
            % Nothing special to do for PECh!
            % Do nothing for PBC either.
        otherwise
            error(['Unknown right boundary condition type ' num2str(BCr)]);
    end
    switch(BCl)                                                             % [2006-11-07] COMPLETED
        case PBC
            % Do nothing - weighting factors for fields below boundary must come from top of domain
            Hzr_ldp(1,:) = Hzr_ld(1,:); Hzr_ld(1,:) = 0;                    % [required] Hzr_ldp and Hzr_lup = new matrix diagonals for left PBC wrapping
            Hzr_lup(1,:) = Hzr_lu(1,:); Hzr_lu(1,:) = 0;                    % [required]
        case {PMC,PMCh}
            Hzr_lu(1,:) = 0; Hzr_ld(1,:) = 0;                               % [required] left row of normal E-fields wrap around to right if not zeroed
            if(BCl==PMC), Hzr_ru(1,:) = 2*Hzr_ru(1,:); Hzr_rd(1,:) = 2*Hzr_rd(1,:); end % [required] for PMC needs doubling, for PMCh stays the same
        case {PEC,PECh}
            Hzr_lu(1,:) = 0; Hzr_ld(1,:) = 0; Hzr_ru(1,:) = 0; Hzr_rd(1,:) = 0; % [safety?] left columns of normal E-fields correspond to an output field that is zero and will be removed from matrix below
            % Do nothing special for PECh.
        otherwise
            error(['Unknown left boundary condition type ' num2str(BCl)]);
    end
    %--- End of setting Hzr boundary conditions

%% Confirmed not needed.  Works fine before adding this:
%    Hzr_lu = [Hzr_lu, zeros(M,1)]; Hzr_ru = [Hzr_ru, zeros(M,1)]; Hzr_lup = [Hzr_lup, zeros(M,1)];          % Fill top with zeros to match the input Ez field in N dimension size before diagonalizing block
%    Hzr_ld = [Hzr_ld, zeros(M,1)]; Hzr_rd = [Hzr_rd, zeros(M,1)]; Hzr_ldp = [Hzr_ldp, zeros(M,1)];

    Hzr = spdiags2([ Hzr_lu(:)*(BCu==PBC) Hzr_ru(:)*(BCu==PBC) Hzr_ld(:) Hzr_rd(:) Hzr_ldp(:)*(BCr==PBC) Hzr_lu(:) Hzr_ru(:) Hzr_lup(:)*(BCr==PBC)], ...  % Complete BC support in Hzr
                  [       (N-2)*M+1            (N-2)*M           1         0              -(M-2)         -(M-1)      -M          -2*(M-1)        ], ...
                        M*N, M*(N-1)).';
% input size is (M-1 x N) -^      ^- output size (M x N-1); note: input size has to be decimated by removing columns
    Hzr((N-2)*M+1,M-1) = Hzr_lup(1,N-1)*(BCr==PBC)*(BCu==PBC); % Single matrix element for top left corner wrapping to bottom right corner when both PBCs are on.
    clear Hzr_ld Hzr_lu Hzr_rd Hzr_ru Hzr_ldp Hzr_lup;

    nidx = setdiff((1:M*N), M*(1:N));
    Hzr = Hzr(:,nidx);      % [MP] Reformat the Hzr matrix eliminating the above empty lines.
%%%end

%    spy(Hzr);
else    % For Semi-Vectorial computation (OPTS.fieldmode == 'M?'), Hrz = Hzr = 0.
    Hrz = sparse(size(Hrr,1),size(Hzz,2)); Hzr = Hrz.';  % Empty matrices
end


% Subtract either (w^2 u0 e0 er) or (beta^2) (i.e. (gamma/rho)^2):         % [2006-11-03] SECTION COMPLETED
if OPTS.eigmode == 'w' % In this case, keig = gamma^2 = (beta rho)^2 = (number of wavelengths around ring)^2
%    A = rho_rr(:,2:N-1).^2; Hrr = Hrr + spdiags(1./A(:),0,size(Hrr,1),size(Hrr,2)) * keig; clear A;
%    A = rho_zz(2:M-1,:).^2; Hzz = Hzz + spdiags(1./A(:),0,size(Hzz,1),size(Hzz,2)) * keig; clear A;
%    A = er_rr(:,2:N-1); A = spdiags(1./A(:),0,size(Hrr,1),size(Hrr,2)); % Multiply by inverse square index!
%    Hrr = A * Hrr; Hrz = A * Hrz; clear A;
%    A = er_zz(2:M-1,:); A = spdiags(1./A(:),0,size(Hzz,1),size(Hzz,2));
%    Hzz = A * Hzz; Hzr = A * Hzr; clear A;
    A = rho_rr(:,1:N).^2; Hrr = Hrr + spdiags2(1./A(:),0,size(Hrr,1),size(Hrr,2)) * keig; clear A;
    A = rho_zz(1:M,:).^2; Hzz = Hzz + spdiags2(1./A(:),0,size(Hzz,1),size(Hzz,2)) * keig; clear A;
    A = er_rr(:,1:N); A = spdiags2(1./A(:),0,size(Hrr,1),size(Hrr,2)); % Multiply by inverse square index!
    Hrr = A * Hrr; Hrz = A * Hrz; clear A;
    A = er_zz(1:M,:); A = spdiags2(1./A(:),0,size(Hzz,1),size(Hzz,2));
    Hzz = A * Hzz; Hzr = A * Hzr; clear A;
    % Eigenvalues of H will be frequencies normalized to speed of light, (w/c)^2
    if (VMODE ~= 0), fprintf('m2dpmloper: w-eigenvalue\n'); end
else    % Assume eigmode == 'b'; then, keig = k0^2 = w^2 u0 e0
%    A = rho_rr(:,2:N-1).^2; B = er_rr(:,2:N-1); NH = size(Hrr);
%    Hrr =  spdiags(A(:),0,NH(1),NH(2)) * (-Hrr + spdiags(B(:),0,NH(1),NH(2)) * keig);
%    Hrz = -spdiags(A(:),0,NH(1),NH(2)) * Hrz; clear A B NH;
%    A = rho_zz(2:M-1,:).^2; B = er_zz(2:M-1,:); NH = size(Hzz);
%    Hzz =  spdiags(A(:),0,NH(1),NH(2)) * (-Hzz + spdiags(B(:),0,NH(1),NH(2)) * keig);
%    Hzr = -spdiags(A(:),0,NH(1),NH(2)) * Hzr; clear A B NH;
%%    Hrz = -Hrz; Hzr = -Hzr;
    A = rho_rr(:,1:N).^2; B = er_rr(:,1:N); NH = size(Hrr);
    Hrr =  spdiags2(A(:),0,NH(1),NH(2)) * (-Hrr + spdiags2(B(:),0,NH(1),NH(2)) * keig);
    Hrz = -spdiags2(A(:),0,NH(1),NH(2)) * Hrz; clear A B NH;
    A = rho_zz(1:M,:).^2; B = er_zz(1:M,:); NH = size(Hzz);
    Hzz =  spdiags2(A(:),0,NH(1),NH(2)) * (-Hzz + spdiags2(B(:),0,NH(1),NH(2)) * keig);
    Hzr = -spdiags2(A(:),0,NH(1),NH(2)) * Hzr; clear A B NH;
    % Eigenvalues of H will be propagation constants, gamma^2 = beta^2*rho^2
    if (VMODE ~= 0), fprintf('m2dpmloper: beta-eigenvalue\n'); end
end

%% Code section removes matrix rows/columns for field components that are  % [2006-11-03] MODIFICATION COMPLETED
%% fixed a priori by choice of BCs (i.e. PMCh):
ix = [];
%if (BCl == PMCh) ix = [ix [  1:(M-1):(M-1)*(N-2)] ]; end                  % [MP-21nov2003] Remove left Hrr values
%if (BCr == PMCh) ix = [ix [M-1:(M-1):(M-1)*(N-2)] ]; end                  % [MP-21nov2003] Remove right Hrr values
if (BCd == PEC || BCd == PECh),                ix = [ix (1:1:M-1)                 ];  end  % [2006-11-03-MP] Remove bottom Hrr values
if (BCu == PEC || BCu == PECh || BCu == PBC),  ix = [ix ((M-1)*(N-1) + (1:1:M-1)) ];  end  % [2006-11-03-MP] Remove top Hrr values
ix = setdiff((1:(M-1)*N),ix);                                              % Keep the complementary set
Hrr = Hrr(ix,ix); Hrz = Hrz(ix,:); Hzr = Hzr(:,ix);
ix = [];
%if (BCd == PMCh) ix = [ix [              1:1:(M-2)] ]; end  % [MP-21nov2003] Remove bottom Hzz values
%if (BCu == PMCh) ix = [ix [(M-2)*(N-2) + 1:1:(M-2)] ]; end  % [MP-21nov2003] Remove top Hzz values
if (BCl == PEC || BCl == PECh),                ix = [ix (1:M:M*(N-1)) ];  end  % [2006-11-03-MP] Remove left Hzz values
if (BCr == PEC || BCr == PECh || BCr == PBC),  ix = [ix (M:M:M*(N-1)) ];  end  % [2006-11-03-MP] Remove right Hzz values
ix = setdiff((1:M*(N-1)),ix);
Hzz = Hzz(ix,ix); Hrz = Hrz(:,ix); Hzr = Hzr(ix,:);

H = [Hrr Hrz; Hzr Hzz];
clear Hrr Hrz Hzr Hzz;
%figure; spy(H);

if (VMODE ~= 0)
    fprintf('m2dpmloper: Elapsed time: %f sec\n', etime(clock,t0));
end
%---------.---------.---------.---------.---------.---------.---------.--------|
% END OF MAIN CODE                                                             |
%---------.---------.---------.---------.---------.---------.---------.--------|


%---------.---------.---------.---------.---------.---------.---------.--------|
% mresolve() - Splits large matrix into 3-4 interleaved smaller matrices.      |
%   A    is (1:2M-1)     x (1:2N-1)
%   A_rr is 2:2:(2M-1)-1 x 1:2:2N-1     (includes edges)
%   A_zr is 1:2:(2M-1)   x 1:2:2N-1
%   A_zz is 1:2:(2M-1)   x 2:2:(2N-1)-1
%   A_zr is 2:2:(2M-1)-1 x 2:2:(2N-1)-1
%---------.---------.---------.---------.---------.---------.---------.--------|
function [A_rr, A_zr, A_zz, A_rz] = mresolve(A)
M = (size(A,1)+1)/2; N = (size(A,2)+1)/2; % Size of *basic* grid.
%[MP]Jun28,02
P = size(A,3);                            % [MP] Jun 28,02 Adding tensor index!
if (P == 3), p = [1 2 3]; else p = [1 1 1]; end  % Is refr index tensor or isotr?
%---end of Jun28,02 added code (plus ,p(?) added in the 3 lines below!)

A_rr = A(2*(1:M-1),2*(1:N)-1,p(1));       % (M-1)xN     Er-like matrix [MP]Jun28
A_zr = A(2*(1:M)-1,2*(1:N)-1,p(3));       % MxN         matrix         [MP]Jun28
A_zz = A(2*(1:M)-1,2*(1:N-1),p(2));       % Mx(N-1)     Ez-like matrix [MP]Jun28
if nargout > 3
    A_rz = A(2*(1:M-1),2*(1:N-1));        % (M-1)x(N-1) matrix
else
    A_rz = [];                            % If not asked for, just define it..
end

function res1 = spdiags2(arg1,arg2,arg3,arg4)
%SPDIAGS2 Sparse matrix formed from diagonals.
%   SPDIAGS2 is the repaired version of SPDIAGS, a native Matlab function.
%
%   [NOTE - 2006-11-07 Milos Popovic]:  For rectangular matrices (m ~= n),
%   A = SPDIAGS(B,d,m,n) shifts diagonals up for d>0 and down for d<0; but
%   this works properly only for m>=n.  If m<n, then you get a different
%   matrix; for m<n, diagonals are moved right for d>0 and left for d<0.
%
% Example:
%
% The following code (3 lines) should produce a 7x7, 7x6 or 6x7 chunk of
% the same matrix, but it does not.  For m<n (m=6, n=7), spdiags functions
% qualitatively differently.  This function fixes this problem, so you get
% subsections of the same matrix.
%
% A=[1:7].' * [1 1 1 1]
% H = full( spdiags(A,[-5 -3 0 1],7,7) )   % Moves off diagonals vertically - ok, consistent with linpack
% H = full( spdiags(A,[-5 -3 0 1],7,6) )   % Moves off diagonals vertically - ok, consistent with linpack
% H = full( spdiags(A,[-5 -3 0 1],6,7) )   % Moves off diagonals horizontally - this is inconsistent
%
%

if (arg3 < arg4)
    res1 = spdiags(arg1,arg2,arg4,arg4);
    res1 = res1(1:arg3,:);
else
    res1 = spdiags(arg1,arg2,arg3,arg4);
end
