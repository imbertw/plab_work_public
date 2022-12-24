% Lecture 13 - 2D propagation of waves
clear all; close all;
tic

%% 0. Globals
u0 = 4e-7*pi;   % H/m
c = 299792458;  % m/s
Zo = u0*c;      % Impedance of free space


%% 1. Problem setup
lam0 = 1.2;
k0 = 2*pi/lam0;
xymin = -5; xymax = 5; dxy = 1/40; %0.010;
wPML = 1; NPMLorder = 3; xI1um = 2*0.125/4*10;  % xI1um is absorption strength of PML


%% 2. Set up computational domain with PML absorbing layers on sides
xRvec = [xymin:dxy/2:(xymax + dxy/1e4)].'; xvec = xRvec;            % Half-dx discretization to get both Ez and Hx/Hy grids in one grid

% Generate complex-coordinate domain
xPMLright = xRvec(end)-wPML;  idx = xRvec > xPMLright;  % Find x points for right PML
xIvec(idx,:) = +xI1um*(xRvec(idx)-xPMLright).^NPMLorder;
xPMLleft  = xRvec(1)+wPML;    idx = xRvec < xPMLleft;   % Find x points for left PML
xIvec(idx,:) = -xI1um*abs(xRvec(idx)-xPMLleft).^NPMLorder;
xvec = xRvec - j*xIvec;

yvec = xvec;
xEz = xvec(3:2:end-2); yEz = yvec(3:2:end-2);       % Define staggered subgrids for Ez, Hx and Hy
xHy = xvec(2:2:end-1); yHx = yvec(2:2:end-1);

figure; plot(xRvec, [real(xvec) imag(xvec)], '-', 'LineWidth', 2); grid on;
title('Complex x domain (y domain is the same)'); xlabel('x_R (\mum)'); xlabel('x_R, x_I (\mum)');
legend('x_R','x_I');
set(gcf,'Color',[1 1 1]);

[X,Y] = meshgrid(xvec,yvec); X = X.'; Y = Y.';      % Get X,Y grids and rotate so columns are x direction.


%% 3. Make refractive index distribution

% Low dimension example:
%[X,Y] = meshgrid([1:6],[20:24])
% (X-3).^2 + (Y-22).^2

%k^2 = omega^2*mu*epsilon = omega^2*mu*epsilon0*n^2 = (k0*n)^2
nmat = ones(size(X)) * 1;                                       % Uniform media, air
% nmat(abs((X-1).^2 + (Y-0.5).^2) < 3.^2) = 3;                  % Draw ring resonator
% nmat(abs((X-1).^2 + (Y-0.5).^2) < 2.5.^2) = 1;
% nmat(abs((X-1).^2 + (Y).^2) < 1.1.^2) = 3;                  % Draw ring resonator
% nmat(abs((X-1).^2 + (Y).^2) < 0.9.^2) = 1;
%nmat( X < 0 ) = 3;                    % Draw flat dielectric interface
%nmat( X < -1 ) = 3;                  % Draw flat dielectric interface
%nmat( X < 2 ) = 3;                   % Draw flat dielectric interface
%nmat(abs(Y) < 0.1) = 3;                  % Draw waveguide
% nmat(abs(Y) < 0.1 & X < 1) = 3;                                 % Draw horizontal bus waveguide
% nmat((abs((X-1).^2 + (Y-1).^2) <= 1.1.^2) & X >= 1 & Y <= 1) = 3;                  % Draw bend
% nmat((abs((X-1).^2 + (Y-1).^2) < 0.9.^2) & X >= 1 & Y <= 1) = 1;                  % Draw bend
% nmat(abs(X-2) < 0.1 & Y >= 1) = 3;                                 % Draw horizontal bus waveguide
nmat(abs(Y) < 0.1) = 3;                  % Draw waveguide
nmat(abs((X-1).^2 + (Y-1.3).^2) < 1.1.^2) = 3;                  % Draw ring resonator
nmat(abs((X-1).^2 + (Y-1.3).^2) < 0.9.^2) = 1;
nmat(abs(Y-2.6) < 0.1) = 3;                  % Draw waveguide


figure; imagesc(real(xvec), real(yvec), abs(nmat.'));
set(gca,'ydir','normal'); axis image; colorbar; colormap(gray);
xlabel('Position x (\mum)'); ylabel('Position y (\mum)'); title('Refractive index distribution (at dxy/2 grid), n(x,y)');

nmat2 = nmat(3:2:end-2,3:2:end-2);  % Take index only where Ez field is.

%nvec = nmat2.'; nvec = nvec(:);      % Unwrap 2D structure into column vector
nvec = nmat2(:);                     % Unwrap 2D structure into column vector
figure; imagesc(real(xEz), real(yEz), abs(nmat2.'));
set(gca,'ydir','normal'); axis image; colorbar; colormap(gray);
xlabel('Position x (\mum)'); ylabel('Position y (\mum)'); title('Refractive index distribution (downsampled to Ez''s dxy grid), n(x,y)');
set(gcf,'Color',[1 1 1]);

% figure; plot(abs(nvec), '-o', 'LineWidth', 2);        % Plot unwrapped refractive index distribution
% xlabel('Unwrapped x dimension'); ylabel('Refractive index');
% ylim( [0.9 3.5]);


%% 4. Make source

% Point source
% Jz = double((X == 0) & (Y == 0));

% % 1D Array
% wd=dxy/10;
% xsrc = [-2.25:dxy:2.25];
% ysrc = zeros(size(xsrc));
% Jz = zeros(size(X));
% for kk = 1:length(xsrc)
%     Jz = Jz + exp(-((X-xsrc(kk))/wd).^2) .* exp(-((Y-ysrc(kk))/wd).^2);
% end
%Copy to end to see far field
%idx = find(yEz==19)
%figure; plot(real(xEz), [real(Ezmat(:,idx)) abs(Ezmat(:,idx)) abs(2e-3*sinc(real(xEz)/4.3))], '-', 'LineWidth', 1.5)

% Gaussian beam
% wdaper = 0.5;
% Jz = exp(-(X/wdaper).^2) .* exp(-(Y/dxy).^2);
wdaper = 0.15;
Jz = exp(-(X/dxy).^2) .* exp(-(Y/wdaper).^2);


% figure`; imagesc(real(xvec), real(yvec), abs(Jz));
% set(gca,'ydir','normal'); axis image; colorbar; colormap(gray);
% xlabel('Position x (\mum)'); ylabel('Position y (\mum)'); title('Current density distribution, J_z(x,y)');

Jz2 = Jz(3:2:end-2,3:2:end-2); Jzvec = Jz2(:);
figure; imagesc(real(xEz), real(yEz), abs(Jz2.'));
set(gca,'ydir','normal'); axis image; colorbar; colormap(gray);
xlabel('Position x (\mum)'); ylabel('Position y (\mum)'); title('Current density distribution, J_z(x,y)');
set(gcf,'Color',[1 1 1]);


%% 5. Make wave equation operator

% Using matrix operators to do differentiation:

% 5a. Create d/dx - first first-derivative in x
NN = (length(xvec)-1)/2+1;                      % Domain width for Ez (here equal to height, since we made a square domain)
A = spdiags(-ones(NN^2,1),-1, NN^2, NN^2); A = A(2:end,:);
B = spdiags( ones(NN^2,1),+1, NN^2, NN^2); B = B(1:end-1,:);

% Filter out rows giving unphysical derivatives
idx = [NN+1:NN:NN^2]-1;                                 % First derivative non-existent points
C = setdiff(1:NN^2-1,idx);                                % Keep only all the other rows of A,B
A = A(C,:); B = B(C,:);
% figure; imagesc(A+B); colorbar; colormap(redbluehilight); axis image;
% xlabel('E_z full-grid input'); ylabel('H_y half-grid output');
% title('First first derivative in x - Ez pattern only (full grid in, half grid out)');

% Make non-uniform 1/dx on full grid
idxmatrixf = 1./(X(3:2:end,1:2:end)-X(1:2:end-2,1:2:end));          % Non-uniform grid 1/dx at all Full grid x points, Full grid y points
idxmatrixf = spdiags(idxmatrixf(:), 0, (NN-1)*NN, (NN-1)*NN);
DXmatf = idxmatrixf * (A+B);                            % First derivative operator
% figure; imagesc(real(DXmatf)); colorbar; colormap(redbluehilight); axis image;
% xlabel('E_z full-grid input'); ylabel('H_y half-grid output');
% title('First first derivative in x - full matrix (full grid in, half grid out)');


% 5b. Second first derivative
A = spdiags(-ones(NN^2,1),-1, (NN-1)*NN, (NN-1)*NN); A = A(2:end,:);     % Now make dx for half-grid
B = spdiags( ones(NN^2,1),+1, (NN-1)*NN, (NN-1)*NN); B = B(1:end-1,:);

% Filter out rows giving unphysical derivatives
idx = [NN+1-1:NN-1:(NN-1)*NN]-1;                        % First derivative non-existent points
C = setdiff(1:(NN-1)*NN-1,idx);                         % Keep only all the other rows of A,B
A = A(C,:); B = B(C,:);
% figure; imagesc(A+B); colorbar; colormap(redbluehilight); axis image;
% xlabel('H_y half-grid input'); ylabel('E_z half-grid output');
% title('Second first derivative in x (half grid in, full grid out)');

% Make non-uniform 1/dx on full grid
idxmatrixh = 1./(X(4:2:end-1,1:2:end)-X(2:2:end-3,1:2:end));        % Non-uniform grid 1/dx at all Half grid x points, Full grid y points
idxmatrixh = spdiags(idxmatrixh(:), 0, (NN-2)*NN, (NN-2)*NN);
DXmath = idxmatrixh * (A+B);      % First derivative operator

DX2mat = DXmath * DXmatf;
% figure; imagesc(real(DX2mat)); colorbar; colormap(redbluehilight); axis image;
% xlabel('E_z full-grid input'); ylabel('E_z full-grid (2 pixels smaller) output');
% title('Second derivative in x (full grid in, full grid out)');

% Drop out unused (zero by default) fields
idx = [1:NN NN^2-NN+1:NN^2 NN+1:NN:NN^2-NN 2*NN:NN:NN^2-NN];
C = setdiff(1:NN^2,idx);                         % Keep only all the other rows of A,B
DX2mat = DX2mat((NN-2)+1:end-(NN-2),C);
% figure; imagesc(real(DX2mat)); colorbar; colormap(redbluehilight); axis image;
% xlabel('E_z full-grid input'); ylabel('E_z full-grid (2 pixels smaller) output');
% title('Second derivative in x (full grid in, full grid out)');


% 5c. Create d/dy - first first-derivative in y
NN = (length(xvec)-1)/2+1;                      % Domain width for Ez (here equal to height, since we made a square domain)
A = spdiags(-ones(NN^2,1), 0,   NN^2, NN^2); A = A(1:end-NN,:);
B = spdiags( ones(NN^2,1), +NN, NN^2, NN^2); B = B(1:end-NN,:);
% figure; imagesc(A+B); colorbar; colormap(redbluehilight); axis image;
% xlabel('E_z full-grid input'); ylabel('H_x half-grid output');
% title('First first derivative in y - Ez pattern only (full grid in, half grid out)');

% Make non-uniform 1/dx on full grid
idymatrixf = 1./(Y(1:2:end,3:2:end)-Y(1:2:end,1:2:end-2));          % Non-uniform grid 1/dy at all Full grid y points
idymatrixf = spdiags(idymatrixf(:), 0, NN*(NN-1), NN*(NN-1));
DYmatf = idymatrixf * (A+B);                            % First derivative operator
% figure; imagesc(real(DYmatf)); colorbar; colormap(redbluehilight); axis image;
% xlabel('E_z full-grid input'); ylabel('H_y half-grid output');
% title('First first derivative in y - full matrix (full grid in, half grid out)');


% 5d. Second first derivative in y
A = spdiags(-ones(NN^2,1), 0,   NN*(NN-1), NN*(NN-1)); A = A(1:end-NN,:);     % Now make dy for half-grid
B = spdiags( ones(NN^2,1), +NN, NN*(NN-1), NN*(NN-1)); B = B(1:end-NN,:);
% figure; imagesc(A+B); colorbar; colormap(redbluehilight); axis image;
% xlabel('H_x half-grid input'); ylabel('E_z half-grid output');
% title('Second first derivative in y (half grid in, full grid out)');

% Make non-uniform 1/dx on full grid
%idymatrixh = 1./(Y(1:2:end,4:2:end-1)-Y(1:2:end,1:2:end-3));        % Non-uniform grid 1/dx at all Half grid x points, Full grid y points
idymatrixh = 1./(Y(1:2:end,4:2:end-1)-Y(1:2:end,2:2:end-3));        % Non-uniform grid 1/dx at all Half grid x points, Full grid y points
idymatrixh = spdiags(idymatrixh(:), 0, NN*(NN-2), NN*(NN-2));
DYmath = idymatrixh * (A+B);      % First derivative operator

DY2mat = DYmath * DYmatf;
% figure; imagesc(real(DY2mat)); colorbar; colormap(redbluehilight); axis image;
% xlabel('E_z full-grid input'); ylabel('E_z full-grid (2 pixels smaller) output');
% title('Second derivative in y (full grid in, full grid out)');

% Drop out unused (zero by default) fields
DY2mat = DY2mat(:,NN+1:end-NN);
idx = [1:NN:NN*(NN-2) NN:NN:NN*(NN-2)];
C = setdiff(1:NN*(NN-2),idx);                         % Keep only all the other rows of A,B

DY2mat = DY2mat(C,C);
% figure; imagesc(real(DY2mat)); colorbar; colormap(redbluehilight); axis image;
% xlabel('E_z full-grid input'); ylabel('E_z full-grid (2 pixels smaller) output');
% title('Second derivative in y (full grid in, full grid out)');

% 5e:  Laplacian in xy
DXY2mat = DX2mat+DY2mat;
% figure; imagesc(real(DXY2mat)); colorbar; colormap(redbluehilight); axis image;
% c = caxis; caxis([-1 1]*max(abs(c(:))));
% xlabel('E_z full-grid input'); ylabel('E_z full-grid (2 pixels smaller) output');
% title('Laplacian in xy (full grid in, full grid out)');


% 6: Make full wave equation operator

Hmat = DXY2mat + spdiags( (k0 * nvec).^2, 0, (NN-2)^2, (NN-2)^2);
% figure; imagesc(abs(Hmat)); colorbar; colormap(redbluehilight); axis image;
% c = caxis; caxis([-1 1]*max(abs(c(:))));
% Hmax0 = max(abs(Hmat(:))); caxis([-1 1]*Hmax0);
disp('Generated all the operators...');

%% Solve for the radiated field:

%Ez = inv(Hmat) * Jz;
Ezvec = Hmat \ Jzvec;
disp('Solved the matrix problem...');
toc

Ezmat = reshape(Ezvec, NN-2, NN-2);

figure; imagesc(real(xEz), real(yEz), real(Ezmat.'));
set(gca,'ydir','normal'); axis image; colorbar; colormap(redbluehilight);
c = caxis; caxis([-1 1]*max(abs(c(:))));
xlabel('Position x (\mum)'); ylabel('Position y (\mum)'); title('Ez field distribution, E_z(x,y)');
set(gcf,'Color',[1 1 1]);

figure; imagesc(real(xEz), real(yEz), abs(Ezmat.'));
set(gca,'ydir','normal'); axis image; colorbar; colormap(redbluehilight);
xlabel('Position x (\mum)'); ylabel('Position y (\mum)'); title('Ez field distribution, E_z(x,y)');
set(gcf,'Color',[1 1 1]);

% Ez movie
wt=pi/2;
%figure;
%plot(xvec, real(Ez*exp(j*wt)), xvec, imag(Ez*exp(j*wt)), 'LineWidth', 2);

Ezmax = max(abs(Ezmat(:)));
figure;
for mm = 1:24
    wt = (mm-1)/24 * 2*pi;
%    plot(xRvec, [real(Ez*exp(j*wt)) abs(Ez) -abs(Ez)], 'LineWidth', 2);
%    grid on; ylim([-1.5 1.5]*Ezmax); xlabel('Position x (\mum)'); ylabel('E_z field (a.u.)');
%     set(gca,'ydir','normal'); axis image; colorbar; colormap(redbluehilight);
    imagesc(real(xEz), real(yEz), real(Ezmat.'*exp(j*wt)));
    set(gca,'ydir','normal'); axis image; colorbar; colormap(redbluehilight);
    caxis([-1 1]*Ezmax);
    xlabel('Position x (\mum)'); ylabel('Position y (\mum)'); title('E_z radiation given J_z source (2D)');
%    xlim(auto);
    set(gcf,'Color',[1 1 1]);
    
    M(mm) = getframe;
end
movie(M,5)
