%3DField.m
function [E_3D,X,Z,ratio]=ThreeDFieldLong(a,L,n1,n2,n3,h1,h2,q1,q2,beta1,beta2,ZE1,ZE2,ZE3,PhiE_mirror,PhiE_cavity,x_x,y_y,z_z)
%Generate  3D longitudinal mode
clear ZE1_3D ZE2_3D ZE3_3D
ZE1_3D(1,1,:)=ZE1; %first, make a 1x1xsomething array with values of the longitudinal mode
ZE1_3D = repmat(ZE1_3D,[size(PhiE_mirror,1) size(PhiE_mirror,2) 1]); %take those values and spread them over a 2D plane the size of the cross-section
ZE2_3D(1,1,:)=ZE2;
ZE2_3D = repmat(ZE2_3D,[size(PhiE_cavity,1) size(PhiE_cavity,2) 1]);
ZE3_3D(1,1,:)=ZE3;
ZE3_3D = repmat(ZE3_3D,[size(PhiE_mirror,1) size(PhiE_mirror,2) 1]);
%Generate 3D transverse mode
PhiE_mirror_3D = repmat(PhiE_mirror,[1 1 length(ZE1)]); %make a 3D array which is the transverse mode repeated in the z direction
PhiE_cavity_3D = repmat(PhiE_cavity,[1 1 length(ZE2)]);
%Multiply to obtain full 3D field
E_3D_1= PhiE_mirror_3D.*ZE1_3D; %multiply transverse and longitudinal values elementwise to produce field
E_3D_2 = PhiE_cavity_3D.*ZE2_3D;
E_3D_3 = PhiE_mirror_3D.*ZE3_3D;
%concatenate
E_3D = cat(3,E_3D_1,E_3D_2,E_3D_3); %concatenate in the z dimension
%Generate coordinate maps
[X,Z] = meshgrid(z_z,x_x);
Elong = E_3D(100,:,:);
Elong = squeeze(Elong);
ratio = a/L;

%Normalize Fields

B = besselj(1,h1*a)/besselk(1,q1*a);
W = cos(beta1*L/2)/exp(1i*beta2*L/2);
intepsE_core_cavity = a*n1^2*pi*(a*h1*besselj(0,a*h1)^2-2*besselj(0,a*h1)*besselj(1,a*h1)+a*h1*besselj(1,a*h1)^2)*(beta1*L+sin(beta1*L))/(2*beta1*h1);
intepsE_clad_cavity = B^2*n2^2*pi^(1.5)*meijerG([],[1.5],[0 0 2],[],a^2*q1^2)*(beta1*L+sin(beta1*L))/(4*beta1*q1^2);
intepsE_core_mirror = 1i*a*exp(1i*beta2*L)*n2^2*pi*W^2*(a*h2*besselj(0,a*h2)^2-2*besselj(0,a*h2)*besselj(1,a*h2)+a*h2*besselj(1,a*h2)^2)/(beta2*h2);
intepsE_clad_mirror = 1i*B^2*exp(1i*beta2*L)*n3^2*pi^(1.5)*W^2*meijerG([],[1.5],[0 0 2],[],a^2*q2^2)/(2*beta2*q2^2);

Elong = Elong/sqrt(intepsE_core_cavity+intepsE_clad_cavity+intepsE_core_mirror+intepsE_clad_mirror);
%--------------------------------------------------------------------------
%Plot Fields
%--------------------------------------------------------------------------
set(gcf,'color','w');
colormap(redhilight)
pcolor(Z.',X.',(Elong).')
caxis([0 1.4])
ylabel('z/\mum')
xlabel('x/\mum')
colorbar
shading interp
hold on
theta = meshgrid(linspace(0, 2*pi, 50), linspace(0, 2*pi, 50)) ;

newX = a .* cos(theta);
newY = a .* sin(theta);

newZ = meshgrid(linspace(-L/2, L/2, 50), linspace(-L/2, L/2, 50))';

surf(newY.',newZ.',newX.','FaceColor','k','FaceAlpha',0.2,'LineStyle','none')

hold off
axis image


