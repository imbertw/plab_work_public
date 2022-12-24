
%This is a simple example of using the modesolver to find the modes of a 
%Single Si waveguide on SiO2, no overcladding


clear;

%Add the modesolver into the path
addpath('C:\Users\undergrad\Code\Modesolver Updated 20120430');
um = 1;
c = 299792458;

%Material parameters
nSi = 3.47;
nAir = 1;
nSiO2 = 1.45;

nlyrs = ones(3, 3);
nlyrs(3,:) = nSiO2;
nlyrs(2,1) = nSiO2; nlyrs(2,3) = nSiO2;
nlyrs(2,2) = nSi;
nlyrs(1,:) = nSiO2;


%Layer thickness and widths
Siwidth = .48;            %waveguide width
Sithick = .220;          %waveguide height
side = 1.5;
top =1.5;
bot = 1.5;
rOut = 12; 

dlyrsx = [side Siwidth side];         %microns wide            
dlyrsy = [bot Sithick top];           %microns tall

%Run modesolver
dxy = [0.01, 0.01];                   %Discretization
k0 = 126;                             %Don't use k0 anymore it is now the modal order
OPTS.mu_guess = 2*pi/(1.55);          %guess of the modal order;
OPTS.NMODES_CALC = 5;                 %Number of modes you want to calculate

OPTS.eigmode = 'w';                   %Omega mode
OPTS.PMLwidth = [0 0.5 0 0];
OPTS.PMLsigma = [0.2 0.2];
OPTS.radius = (rOut-side-Siwidth); 



%The matlab function sisolver3d calls the modesolver
%  [N,F, V] = sisolver3d2(nlyrs, dlyrsx, dlyrsy, dxy, k0, OPTS);
%It outputs three structures
%     N
%     F
%     V  (don't worry about this one)

%N is a struct that contains information about the structure and material
%  N.x = 0:dx/2:width of computational space                                 If you are interested the dx/2 is used instead of dx due to the Yee formulation of discrete EM theory (see Chew paper)
%  N.y = 0:dy/2:height of computational space
%  N.n = matrix of material index 

%So if you want to look at your structure before calculating the modes just do:
N = sisolver3d2(nlyrs, dlyrsx, dlyrsy, dxy, k0); 
imagesc(N.y, N.x, (N.n).'); set(gca, 'YDir', 'normal'); colorbar;


%F is a structure that contains alot of information outputted by the modesolver
%   The most useful are:
%       The mode fields stored in the matrices: 
%           F.Ex(:,:,1), F.Ey(:,:,1), F.Ez(:,:,1), F.Hx(:,:,1), F.Hy(:,:,1), F.Hz(:,:,1) for the first mode
%           F.Ex(:,:,2), F.Ey(:,:,2), F.Ez(:,:,2), F.Hx(:,:,2), F.Hy(:,:,2), F.Hz(:,:,2) for the second mode
%             etc... note that when I say first mode I don't mean fundamental, simply the first mode calculated by the modesolver
%       The propagation constants stored in: F.beta:
%           F.beta(1) is the prop const of the first mode
%           F.beta(2) for the second 
%             etc...

%Actually call the modesolver and find the modes
[N,F] = sisolver3d2(nlyrs, dlyrsx, dlyrsy, dxy, k0, OPTS);



SZ.F = F; SZ.N = N; Neff = F.beta/k0;   %Just stores output of modesolver into a structure that can be read by modeview


%Program to look at the modes
modeview(SZ);



 

 