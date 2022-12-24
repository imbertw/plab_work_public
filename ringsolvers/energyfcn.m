% solver for the total energy of a field profile in a ring cavity
% Savvy, 06/02/2011
% input: N and F are refractive index matrix and field profile from the
%           output of modesolver. 
%        nmode:  which mode to solve
% output: energy, the total energy of the field profile
% implicit requirement: PEC boundary conditions were used; 

function [energy] = energyfcn(N, F, nmode)

% if isempty(nmode) nmode=1;  end

% epsir_x = ones((size(N.n,1)-1)/2,(size(N.n,2)-3)/2 );    %relative dielectric constant
% epsir_y = ones((size(N.n,1)-3)/2,(size(N.n,2)-1)/2 );
epsir_x = N.n(2:2:end-1, 3:2:end-2).^2; 
epsir_y = N.n(3:2:end-2, 2:2:end-1).^2; 
epsir_z = N.n(3:2:end-2, 3:2:end-2).^2; 

epsilon_0 = 8.85e-12; 
mu_0 = 4*pi*1e-7; 

% total EM energy in straight waveguide: 
% % electric energy
% energy1 = 1/4* epsilon_0*(sum(sum(epsir_x.*abs(F.Ex(:,:,nmode)).^2))+sum(sum(epsir_y.*abs(F.Ey(:,:,nmode)).^2))+sum(sum(epsir_z.*abs(F.Ez(:,:,nmode)).^2)))*F.dx*F.dy; 
% % magnetic energy
% energy2 = 1/4* mu_0*(sum(sum(abs(F.Hx(:,:,nmode)).^2))+sum(sum(abs(F.Hy(:,:,nmode)).^2))+sum(sum(abs(F.Hz(:,:,nmode)).^2)))*F.dx*F.dy; 

% total EM energy in ring cavity: 
Rexmtx = F.Rr(:)*ones(1, size(F.Ex,2)); 
Rhymtx = F.Rr(:)*ones(1, size(F.Hy,2)); 
Rhzmtx = F.Rr(:)*ones(1, size(F.Hz,2)); 
Rhxmtx = F.Rz(:)*ones(1, size(F.Hx,2)); 
Reymtx = F.Rz(:)*ones(1, size(F.Ey,2)); 
Rezmtx = F.Rz(:)*ones(1, size(F.Ez,2)); 

% electric energy
energy1 = 2*pi* 1/4* epsilon_0*(sum(sum(Rexmtx .* epsir_x.*abs(F.Ex(:,:,nmode)).^2))+sum(sum(Reymtx .* epsir_y.*abs(F.Ey(:,:,nmode)).^2))+sum(sum(Rezmtx .* epsir_z.*abs(F.Ez(:,:,nmode)).^2)))*F.dx*F.dy;
% magnetic energy
energy2 = 2*pi* 1/4* mu_0*(intint(Rhxmtx,F.Hx(:,:,nmode))+intint(Rhymtx,F.Hy(:,:,nmode))+intint(Rhzmtx,F.Hz(:,:,nmode)))*F.dx*F.dy;
energy = energy1 +energy2;  

function DoubleIntegral = intint(R, field)
DoubleIntegral = sum(sum(R .* abs(field).^2));
if isnan(DoubleIntegral)
    DoubleIntegral = 0;
end
end
end