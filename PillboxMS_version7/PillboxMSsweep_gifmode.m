%PillboxMS in 3D
%dimensions and parameters
set(0,'DefaultFigureWindowStyle','docked');
%a =0.25;
a = linspace(0.19,0.57,100); %microns
lambda = 1.55; %microns
%Indices
n1 = sqrt(12.085); %Silicon unitless
n2 = sqrt(2.085);  %SiO2 unitless
n3 = PillboxIndex(n1,n2); %Mythical Metal
h = figure;
axis tight manual% this ensures that getframe() returns a consistent size
filename = 'longitudinalXS_radius.gif';
for i=1:length(a)
%Find transverse E modes
%n1 n2
[beta1,neff1,h1,q1,A1,B1,PhiE_cavity,Hr_cavity,Hz_cavity,r_r]=PillboxTV(a(i), lambda, n1,n2); %run PillboxTV function to find transverse fields for n1/n2 layer
%n2 n3
[beta2,neff2,h2,q2,A2,B2,PhiE_mirror,Hr_mirror,Hz_mirror,r_r]=PillboxTV(a(i), lambda, n2,n3); %run PillboxTV function to find transverse fields for n2/n3 layer
%Generate longitudinal mode
[L(i), W(i),ZE1,ZE2,ZE3,z_z] = PillBoxLT(0,beta1,beta2); %run PillBoxLT function to get longitudinal mode
%Generate 3D fields

[Elong,X,Z,ratio(i)]=ThreeDFieldLong(a(i),L(i),n1,n2,n3,h1,h2,q1,q2,beta1,beta2,ZE1,ZE2,ZE3,PhiE_mirror,PhiE_cavity,r_r,r_r,z_z);
drawnow
frame =getframe(h);
im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
%Calculate mode volumes
end