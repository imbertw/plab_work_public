%PillboxMS in 3D
%dimensions and parameters
set(0,'DefaultFigureWindowStyle','docked');
a = linspace(0.1420,0.43,100); %microns
lambda = 1.18; %microns
%Indices
%n1 = sqrt(12.085); %Silicon unitless
%n2 = sqrt(2.085);  %SiO2 unitless
n1 = 3.5240; %Silicon
n2 = 1.4483;  %SiO2
n3 = n2;
n4 = PillboxIndex(n1,n2); %Mythical Metal
for i=1:length(a)
%Find transverse E modes
%n1 n2
[beta1,neff1,h1,q1,A1,B1,PhiE_cavity,Hr_cavity,Hz_cavity,r_r]=PillboxTV(a(i), lambda, n1,n2); %run PillboxTV function to find transverse fields for n1/n2 layer
%n2 n3
[beta2,neff2,h2,q2,A2,B2,PhiE_mirror,Hr_mirror,Hz_mirror,r_r]=PillboxTV(a(i), lambda, n3,n4); %run PillboxTV function to find transverse fields for n2/n3 layer
%Generate longitudinal mode
[L(i), W(i),ZE1,ZE2,ZE3,z_z] = PillBoxLT(0,beta1,beta2); %run PillBoxLT function to get longitudinal mode
%Generate 3D fields
[E_3D,X,Y,Z]=ThreeDField(a(i),L(i),ZE1,ZE2,ZE3,PhiE_mirror,PhiE_cavity,r_r,r_r,z_z);
%Calculate mode volumes
[mdvolume1(i)] = PillboxMV(a(i),L(i),n1,n2,n3,n4,h1,h2,q1,q2,beta1,beta2);
end
clear PhiE_cavity Hr_cavity Hz_cavity r_r PhiE_mirror Hr_mirror Hz_mirror W ZE1 ZE2 ZE3 z_z E_3D X Y Z
a2 = linspace(0.1878,0.57,100); %microns
lambda = 1.55; %microns
%Indices
n1 = sqrt(12.085); %Silicon unitless
n2 = sqrt(2.085);  %SiO2 unitless
n3 = n2;
n4 = PillboxIndex(n1,n2); %Mythical Metal
for i=1:length(a2)
%Find transverse E modes
%n1 n2
[beta1,neff1,h1,q1,A1,B1,PhiE_cavity,Hr_cavity,Hz_cavity,r_r]=PillboxTV(a2(i), lambda, n1,n2); %run PillboxTV function to find transverse fields for n1/n2 layer
%n2 n3
[beta2,neff2,h2,q2,A2,B2,PhiE_mirror,Hr_mirror,Hz_mirror,r_r]=PillboxTV(a2(i), lambda, n3,n4); %run PillboxTV function to find transverse fields for n2/n3 layer
%Generate longitudinal mode
[L2(i), W(i),ZE1,ZE2,ZE3,z_z] = PillBoxLT(0,beta1,beta2); %run PillBoxLT function to get longitudinal mode
%Generate 3D fields
[E_3D,X,Y,Z]=ThreeDField_noplot(a2(i),L2(i),ZE1,ZE2,ZE3,PhiE_mirror,PhiE_cavity,r_r,r_r,z_z);
%Calculate mode volumes
[mdvolume2(i)] = PillboxMV(a2(i),L2(i),n1,n2,n3,n4,h1,h2,q1,q2,beta1,beta2);
end

figure
plot([a],[mdvolume1],'LineWidth',2)
hold on
plot([a2],[mdvolume2],'LineWidth',2)
grid on
ylabel('Mode volume(\mum^3)')
hold on
xlim([0.1420 0.42])
ylim([0 0.16])
plot([0.1878 0.1878],[0 0.16],'k')
hold on
plot(0.2514,0.02405,'ko')
hold on
plot(0.1868,0.01019,'ko')
yyaxis right
plot(a,L,'LineWidth',2)
hold on
plot(a2,L2,'LineWidth',2)
hold on
ylabel('Thickness of core(\mum)')
xlabel('Radius of core(\mum)')
set(gcf,'color','w');
set(gca,'FontSize',24);
set(gca,'FontName','Arial');
