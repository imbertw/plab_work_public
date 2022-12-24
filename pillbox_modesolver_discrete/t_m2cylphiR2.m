%s_m2cylphiR2.m: Discrete cylindrical modesolver for the E_phi field in 2D
% Imbert Wang 11/27/2022
tic
%% 0. Globals
%-------------------------------------------------------------------------%
% Global constants to be used in the code
%-------------------------------------------------------------------------%
%micro = 1e-6;   % units for microns in meters
%% 1. Set up 2D index structure
%-------------------------------------------------------------------------%
% Use a similar method to sisolver3d to define the index structure
%-------------------------------------------------------------------------%
% A. Define refractive indices
n1 = sqrt(12.085);            % Index of silicon at 1550nm
n2 = sqrt(2.085);             % Index of silicon dioxide at 1550nm
n3 = PillboxIndex(n1,n2);     % Index of mythical metal that satisfies permittivity condition

% B. Define index regions: layers of refractive indices
nlyrs = ones(3, 2);
nlyrs(3,1) = n2; nlyrs(3,2) = n3;  % Top layer
nlyrs(2,1) = n1; nlyrs(2,2) = n2;  % Middle layer
nlyrs(1,1) = n2; nlyrs(1,2) = n3;  % Bottom layer

% C. Layer thicknesses and radii
radius_core  = 0.2514;         % Pillbox core radius
thick_core  = 0.19193;         % Pillbox core thickness
side = 1.0;                    % Extent of sides
topClad  = 1.0;                % Top Pillbox cladding thickness
botClad = 1.0;                 % Bottom Pillbox cladding thickness

% D. Establishing the layer thicknesses and widths
dlyrsrho = [radius_core  side];         % microns in radius (i.e., half width)            
dlyrsz = [botClad thick_core topClad];           % microns tall

%% 2. Set up coordinate system
%-------------------------------------------------------------------------%
% Before we proceed, we must define the coordinate system
%-------------------------------------------------------------------------%
% A. Radial coordinates

rho_max_init = sum(dlyrsrho);       % Maximum radius, initially--we will round to ensure whole number of points
x_max_init = 0.5*rho_max_init^2;    % Maximum x, initially--we will round to ensure whole number of points
N_rho_core = 100;                   % Number of points in the radial direction for the core of the pillbox
x_core = 0.5*radius_core^2;         % "Radius" of core in x coordinate system
dx = x_core/N_rho_core;             % Discretization for the x-coordinate
N_rho_total = round(x_max_init/dx); % Find a whole number of points for the entire radial coordinates
x_max = N_rho_total*dx;             % The actual maximum x-coordinate
xvec = 0:dx:x_max;                  % x-coordinate, where x = 0.5*rho^2 in 1D
rhovec = (2*xvec).^(1/2);           % Radial coordinate in 1D
rho_max = rhovec(end);              % The actual maximum radius for the simulation

% B. Z- coordinates

z_max_init = sum(dlyrsz);           % Maximum z, initially--we will round to ensure whole number of points
N_z_core = 100;                     % Number of points in the z-direction for the core of the pillbox
z_core = thick_core;                % Thickness of the core
dz = z_core/N_z_core;               % Discretization in the z-direction
N_z_total = round(z_max_init/dz);   % Find a whole number of points for the entire z-coordinate system
z_max = N_z_total*dz;               % The actual maximum z-coordinate
zvec = 0:dz:z_max;                  % z-coordinate in 1D

N_z_total = N_z_total + 1;          % Update number of points, because we start from 0 instead of dz
N_rho_total = N_rho_total + 1;      % Update number of points, because we start from 0 instead of dx

% C. rho-z coordinate grid
[x_x,z_z] = meshgrid(xvec,zvec); % x- and z- coordinates in 2D
rho_rho = (2*x_x).^(1/2);

%% 3. Generate 2D index structure
%-------------------------------------------------------------------------%
% Having defined the index structure and its coordinate system, we can now
% generate the 2D index structure
%-------------------------------------------------------------------------%
n_x_z = ones(N_z_total, N_rho_total);     % Initialize index matrix
dlyrsx = 0.5*dlyrsrho.^2;
for ii=1:length(dlyrsx)
    for jj=1:length(dlyrsz)
        if (ii == 1) && (jj == 1)
            n_x_z(zvec<=dlyrsz(jj), xvec<=dlyrsx(ii)) = nlyrs(jj,ii);% At the 1,1 position, need to make sure that you don't access -1 variable
        elseif (ii==length(dlyrsx)) && (jj == 1)
            n_x_z(zvec<=dlyrsz(jj), logical(xvec>sum(dlyrsx(1:ii-1)) & xvec<=xvec(end))) = nlyrs(jj,ii);
        elseif ii==1
            n_x_z(logical(zvec>sum(dlyrsz(1:jj-1)) & zvec<=sum(dlyrsz(1:jj))), xvec<=dlyrsx(ii)) = nlyrs(jj,ii);
        elseif jj==1
            n_x_z(zvec<=dlyrsz(jj), logical(xvec>sum(dlyrsx(1:ii-1)) & xvec<=sum(dlyrsx(1:ii)))) = nlyrs(jj,ii);
        elseif (ii == length(dlyrsx)) && (jj == length(dlyrsz))
            n_x_z(logical(zvec>sum(dlyrsz(1:jj-1)) & zvec<=zvec(end)), logical(xvec>sum(dlyrsx(1:ii-1)) & xvec<=xvec(end))) = nlyrs(jj,ii);  % If at the end of both rho and z, don't be limited by the initial limits, as those have been subsequently modified
        elseif ii == length(dlyrsx)
            n_x_z(logical(zvec>sum(dlyrsz(1:jj-1)) & zvec<=sum(dlyrsz(1:jj))), logical(xvec>sum(dlyrsx(1:ii-1)) & xvec<=xvec(end))) = nlyrs(jj,ii); % If at the end of z, don't be limited by the initial limits, as those have been subsequently modified
        elseif jj == length(dlyrsz)
            n_x_z(logical(zvec>sum(dlyrsz(1:jj-1)) & zvec<=zvec(end)), logical(xvec>sum(dlyrsx(1:ii-1)) & xvec<=sum(dlyrsx(1:ii)))) = nlyrs(jj,ii); % If at the end of rho, don't be limited by the initial limits, as those have been subsequently modified
        else
            n_x_z(logical(zvec>sum(dlyrsz(1:jj-1)) & zvec<=sum(dlyrsz(1:jj))), logical(xvec>sum(dlyrsx(1:ii-1)) & xvec<=sum(dlyrsx(1:ii)))) = nlyrs(jj,ii);% Default, we just add the nlyrs to the index structure
        end
    end
end

%% 4. Create Operators
%-------------------------------------------------------------------------%
% Now let's create the operators for the wave equation:

%[ -2 (1/n) sqrt(x) dx2 sqrt(x) (1/n) - (1/n) dz2 (1/n) ](n E_phi) = k_0^2 (n E_phi)

% Derived in Google Drive>Shared Drives>Popovic
% Lab>04_Users>Imbert>202111_pillbox_modesolver.pptx
%-------------------------------------------------------------------------%
% A. Create the sparse matrix that will be a template for the other matrices:
N_elements = N_z_total*N_rho_total;
Hsqrtx = sparse(N_elements, N_elements); % Initialize array
% B. Create the sqrt(x) operator

hsqrtx = sqrt(x_x);             % Create square root x operator
hsqrtx = hsqrtx.';              % Transpose and turn into column vector so it can be turned into square matrix 
hsqrtx = hsqrtx(:);             % Transpose and turn into column vector so it can be turned into square matrix
Hsqrtx = spdiags(hsqrtx, 0, Hsqrtx); % Place elements along main diagonal

% C. Create the 1/n operator
N_x_z = sparse(N_elements, N_elements); % Initialize array
Identity = speye(N_elements);      % Identity matrix
n_x_z = n_x_z.';
n_x_z = n_x_z(:);                % Transpose and turn into column vector
N_x_z = spdiags(n_x_z,0, N_x_z); % Place elements along main diagonal
H_invn = Identity/N_x_z;         % Make inverse n array

% D. Create the dx2 operator

% Create tridiagonal matrix
Hdx2 = sparse(N_elements,N_elements);
d = ones(N_elements,1); 
d_pec = d;
d_pec(N_rho_total:N_rho_total:end) = 0; % Enforce PEC boundary conditions. Every N_rho_totalth element is zero, meaning the boundary is zero on the sides
Hdx2 = spdiags(d_pec, -1, Hdx2);
Hdx2 = spdiags(-2*d, 0, Hdx2);
Hdx2 = spdiags(d_pec, +1, Hdx2);
Hdx2 = Hdx2/dx^2;

% E. Create the dz2 operator

Hdz2 = sparse(N_elements, N_elements);
Hdz2 = spdiags(d, -N_rho_total,Hdz2);
Hdz2 = spdiags(-2*d, 0, Hdz2);
Hdz2 = spdiags(d, +N_rho_total, Hdz2);
Hdz2 = Hdz2/dz^2;

%% 5. Create Combined Operator
%-------------------------------------------------------------------------%
% Now let's create the full operator for the wave equation:

%[ -2 (1/n) sqrt(x) dx2 sqrt(x) (1/n) - (1/n) dz2 (1/n) ](n E_phi) = k_0^2 (n E_phi)
%-------------------------------------------------------------------------%
Htot = -2*H_invn*Hsqrtx*Hdx2*Hsqrtx*H_invn - H_invn*Hdz2*H_invn; % Use operators constructed previously to construct total operator

%% 6. Compute the fields and eigenvalues
%-------------------------------------------------------------------------%
% Use matlab's eigenvector/eigenvalue solver to calculate the resonant k0,
% and by extension the resonant wavelength, for this structure.
%-------------------------------------------------------------------------%
k0_guess = 2*pi/1.55;
[V,D] = eigs(Htot,1,k0_guess^2);

k0 = sqrt(diag(D));
lambda_resonant = 2*pi./k0; % Calculate resonant wavelengths of pillbox 

E_phi = N_x_z\V;
E_phi = reshape(E_phi, N_rho_total,N_z_total);
E_phi = E_phi.'; % Calculate fields of pillbox

toc

E_phi_display = [fliplr(E_phi) E_phi];
rho_rho_display = [-fliplr(rho_rho) rho_rho];
z_z_display = [fliplr(z_z) z_z];


% Use interpolation 
rhovec_interp = linspace(min(rhovec),max(rhovec),3000) ; 
[Ri,Zi] = meshgrid(rhovec_interp,zvec) ;
Z1 = griddata(x,y,z1,Ri,Zi) ; 
pcolor(X,Y,Z1)
shading interp

figure;
h = pcolor(rho_rho_display, z_z_display, abs(E_phi_display));
set(h, 'EdgeColor', 'none')
xlabel('\rho/\mum')
ylabel('z/\mum')
colormap("hot")
colorbar('vertical')

n_x_z_display = reshape(n_x_z, N_rho_total, N_z_total);
n_x_z_display = n_x_z_display.';
n_x_z_display = [fliplr(n_x_z_display) n_x_z_display];
figure;
p = pcolor(rho_rho_display, z_z_display,real(n_x_z_display));
set(p,'EdgeColor','none');

% e_phi = V(:,1)./nvec;
% e_phi = e_phi.'*sqrt(u0);
% lambda_vec = 2*pi./sqrt(diag(D));



