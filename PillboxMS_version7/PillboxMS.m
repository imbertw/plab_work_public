%PillboxMS in 3D
%Dimensions and parameters
a = 0.2514; %microns
%1550nm
lambda = 1.55; %microns
%Indices
n1 = sqrt(12.085); %Silicon unitless
n2 = sqrt(2.085);  %SiO2 unitless
n3= n2;
n4 = PillboxIndex(n1,n2); %Mythical Metal

%{
%1180nm
lambda = 1.18; %microns
%Indices
n1 = 3.5240; %Silicon
n2 = 1.4483;  %SiO2
n3= n2;
n4 = PillboxIndex(n1,n2); %Mythical Metal
%}
%%
%Find transverse E modes
%n1 n2: cavity
[beta1,neff1,h1,q1,A1,B1,PhiE_cavity,Hr_cavity,Hz_cavity,r_r]=PillboxTV(a, lambda, n1,n2); %run PillboxTV function to find transverse fields for n1/n2 layer
%n2 n3: mirror
[beta2,neff2,h2,q2,xx,xx,PhiE_mirror,Hr_mirror,Hz_mirror,r_r]=PillboxTV(a, lambda, n3,n4); %run PillboxTV function to find transverse fields for n2/n3 layer
clear xx
%Generate longitudinal mode
[L,W,ZE1,ZE2,ZE3,z_z] = PillBoxLT(0,beta1,beta2); %run PillBoxLT function to get longitudinal mode
%Generate longitudinal index structure
n_z = zeros(1,length(z_z));
n_z(abs(z_z)<L/2) = max([ZE1 ZE2 ZE3]);
n_z(abs(z_z)>=L/2) = min([ZE1 ZE2 ZE3]);
%Plot longitudinal mode
figure
plot(z_z,[ZE1 ZE2 ZE3],'b','LineWidth',2)
hold on
plot(z_z,n_z,'r','LineWidth',1.5)
title('Longitudinal confinement')
xlabel('z/\mum')
ylabel('Field/[AU]')
ylim([min([ZE1 ZE2 ZE3]) 1.1*max([ZE1 ZE2 ZE3])])
%%
% %Plot radial field
% figure
% imagesc(r_r,r_r,PhiE_cavity.^2)
% colormap(hot)
% colorbar
% axis equal
% ylabel('y/\mum')
% xlabel('x/\mum')
% set(gcf,'color','w');
% set(gca,'FontSize',12);
% set(gca,'FontName','Arial');
%Plot radial field
figure
plot(r_r,PhiE_cavity(:,ceil(end/2)).^2)
axis equal
ylabel('y/\mum')
xlabel('x/\mum')
set(gcf,'color','w');
set(gca,'FontSize',12);
set(gca,'FontName','Arial');
%Generate 3D fields
%[E_3D,X,Y,Z]=ThreeDField(a,L,ZE1,ZE2,ZE3,PhiE_mirror,PhiE_cavity,r_r,r_r,z_z);
%Generate 3D fields sweeping along z (for presentation gif)
[E_3D,X,Y,Z]=ThreeDFieldSweep(a,L,ZE1,ZE2,ZE3,PhiE_mirror,PhiE_cavity,r_r,r_r,z_z);

%%
%Calculate mode volumes
[mdvolume] = PillboxMV(a,L,n1,n2,n3,n4,h1,h2,q1,q2,beta1,beta2);















