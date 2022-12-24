%radiusvlength
set(0,'DefaultFigureWindowStyle','docked');
a = linspace(0.19,0.57,100); %microns
lambda = 1.55; %microns
%Indices
n1 = sqrt(12.085); %Silicon unitless
n2 = sqrt(2.085);  %SiO2 unitless
n3 = PillboxIndex(n1,n2); %Mythical Metal

for i=1:length(a)
[h,q,beta1(i)] = CylindricalParameter(a(i),n1,n2,lambda); %beta in um^-1
[h,q,beta2(i)] = CylindricalParameter(a(i),n2,n3,lambda); %beta in um^-1
[L(i),W,ZE1,ZE2,ZE3,z_z] = PillBoxLT(0,beta1(i),beta2(i));
end
plot(L,a,'LineWidth',2)
xlabel('Length/\mum')
ylabel('Radius/\mum')