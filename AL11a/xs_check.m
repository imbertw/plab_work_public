%xs_check.m: shows the cross section from a xs function.
function xs_check(xs, width, lambda)

[nlyrs, dlyrsx, dlyrsy, left_to_rout] = eval([xs '(width,lambda)']);
dxy           = [0.02, 0.015];                   %Discretization
% Parse the input structure
%[MM,NN] = size(nlyrs);                                                    % Step-index regions
dlyrsx = dlyrsx(:).'; dlyrsy = dlyrsy(:).';
MM = length(dlyrsx); NN = length(dlyrsy);
dx = dxy(1); dy = dxy(2);

dlyrsx(end) = round(sum(dlyrsx)/dx)*dx + 0*1e-3*dx - sum(dlyrsx(1:end-1));
dlyrsy(end) = round(sum(dlyrsy)/dy)*dy + 0*1e-3*dy - sum(dlyrsy(1:end-1)); % [MP-TODO] Here: round xint(end), yint(end) domain *endpoints only* to be even in grid spacing or else round discretization
xint = [0 cumsum(dlyrsx)]; yint = [0 cumsum(dlyrsy)];       % Interface coordinates

% Rounding the coordinates
fprintf('Domain limits: (x,y) = user passed (0..%f, 0..%f), rounding last layer to pixel grid (0..%f, 0..%f) | Pixel grid size = %d x %d\n', sum(dlyrsx), sum(dlyrsy), xint(end), yint(end), round([xint(end)/dx yint(end)/dy]));
fprintf('   x-interface positions at: ');   for k = 1:length(xint), fprintf('%f  ', xint(k)); end;  fprintf('\n');
fprintf('   y-interface positions at: ');   for k = 1:length(yint), fprintf('%f  ', yint(k)); end;  fprintf('\n');

%Grid coordinates
x = (0 : dx/2 : xint(end)); y = (0 : dy/2 : yint(end));     % Grid coordinates
xint(1) = -inf; xint(end) = +inf; yint(1) = -inf; yint(end) = +inf;     % [MP] Put edge interfaces to infinity to avoid index averaging at domain edges

% Standard (area-arithmetic) index averaging for pixels, then those are averaged for pixel edge tensor components for Ex, Ey, Ez
fprintf('Using simple arithmetic dielectric averaging.\n'); 
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

if(find(~isfinite(N.n))), error('sisolver3d error: Generated index distribution has non-finite elements.  Report problem to milos@mit.edu.'); end

%make imaginary values of N.y, N.x and N.n have color 0
Ny = N.y;
Nx = N.x;
Nn = N.n;
Ny(imag(Ny) ~= 0) = 0;
Nx(imag(Nx) ~= 0) = 0;
Nn(imag(Nn) ~= 0) = 0;
imagesc(Ny, Nx, (Nn).'); set(gca, 'YDir', 'normal'); colorbar; % Take a look at the structure