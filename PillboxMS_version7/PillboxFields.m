function [E_phi,E_2D,Hr_2D,Hz_2D,A,B]=PillboxFields(h,q,beta,lambda,a,r1,r2,R)
%--------------------------------------------------------------------------
%Computes fields given cylindrical parameters
%--------------------------------------------------------------------------
%Constants
c = 299792458*1e6;
omega = 2*pi*c/lambda;
mu = 4*pi*1e-13;
%Set normalization constants
A = 1;
B = besselj(1,h*a)/besselk(1,q*a);
%Core Fields (r<a)
E_phi_core = A*besselj(1,h*r1);
H_r_core = -beta/(omega*mu)*A*besselj(1,h*r1);
H_z_core = h*A*r1.*besselj(0,h*r1);
%Cladding Fields (r>a)
E_phi_clad = B*besselk(1,q*r2);
H_r_clad = -beta/(omega*mu)*B*besselk(1,q*r2);
H_z_clad = -q*B*r2.*besselk(0,q*r2);
%Entire field
E_phi1 = [E_phi_core E_phi_clad];
E_phi2 = [fliplr(-E_phi_clad) fliplr(-E_phi_core)];
E_phi = [E_phi2 E_phi1];
%set arrays
E_2D = zeros(size(R));
Hr_2D = zeros(size(R));
Hz_2D = zeros(size(R));
%E_phi field
E_2D(R<=a) = A*besselj(1,h*R(R<=a));
E_2D(R>a) = B*besselk(1,q*R(R>=a));
%H_r field
Hr_2D(R<=a) = -beta/(omega*mu)*A*besselj(1,h*R(R<=a));
Hr_2D(R>a) = -beta/(omega*mu)*B*besselk(1,q*R(R>a));
%H_z field
Hz_2D(R<=a) = h*A*R(R<=a).*besselj(0,h*R(R<=a));
Hz_2D(R>a) = -q*B*R(R>a).*besselk(0,q*R(R>a));
end