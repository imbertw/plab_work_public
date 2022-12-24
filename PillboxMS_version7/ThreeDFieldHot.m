%3DField.m
function [E_3D,X,Y,Z]=ThreeDField(a,L,ZE1,ZE2,ZE3,PhiE_mirror,PhiE_cavity,x_x,y_y,z_z)
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
[X,Y,Z] = meshgrid(x_x,y_y,z_z);

%--------------------------------------------------------------------------
%Plot Fields
%--------------------------------------------------------------------------

figure
colormap(hot)
sx = 0;
sy =0;
sz = z_z(length(z_z)/2);
VolumeSlice = slice(X,Y,Z,abs(E_3D),sx,sy,sz);
%direction = [0 1 0];
%rotate(VolumeSlice,direction,90)
colorbar('vertical');
axis equal
zlabel('z/\mum')
xlabel('x/\mum')
set(gcf,'color','w');
set(gca,'FontSize',24);
set(gca,'FontName','Arial');

%{
for ii = 1:size(VolumeSlice,1);
VolumeSlice(ii,1).FaceColor = 'interp';
VolumeSlice(ii,1).EdgeColor = 'none';
end
%}
VolumeSlice(1,1).FaceColor = 'none';
VolumeSlice(1,1).EdgeColor = 'none';
VolumeSlice(2,1).FaceColor = 'interp';
VolumeSlice(2,1).EdgeColor = 'none';
VolumeSlice(3,1).FaceColor = 'none';
VolumeSlice(3,1).EdgeColor = 'none';

hold on
theta = meshgrid(linspace(0, 2*pi, 50), linspace(0, 2*pi, 50)) ;

newX = a .* cos(theta);
newY = a .* sin(theta);

newZ = meshgrid(linspace(-L/2, L/2, 50), linspace(-L/2, L/2, 50))';

surf(newX,newY,newZ,'FaceColor','k','FaceAlpha',0.2,'LineStyle','none')

view([0 0])

figure
colormap(hot)

%data = E_3D(round(size(E_3D,1)/2),:,:);
%data = squeeze(data);
%data = data.';
%assignin('base','data',data)
%// Define integer grid of coordinates for the above data
%[X,Y] = meshgrid(1:size(data,2), 1:size(data,1));

%// Define a finer grid of points
%[X2,Y2] = meshgrid(1:0.01:size(data,2), 1:0.01:size(data,1));

%// Interpolate the data and show the output
%outData2 = interp2(X, Y, data, X2, Y2, 'cubic');
%imagesc(outData2);

sx = 0;
sy =0;
sz = z_z(length(z_z)/2);
VolumeSlice = slice(X,Y,Z,abs(E_3D),sx,sy,sz);

%direction = [0 1 0];
%rotate(VolumeSlice,direction,90)
%// Cosmetic changes for the axes

% set(gca, 'XTick', linspace(1,size(X2,2),11)); 
% set(gca, 'YTick', linspace(1,size(X2,1),11));
% set(gca, 'XTickLabel', linspace(-0.5,0.5,11));
% set(gca, 'YTickLabel', linspace(-0.5,0.5,11));

%// Add colour bar

colorbar('vertical');
% daspect([2 1 1])
zlabel('z/\mum')
xlabel('x/\mum')
set(gcf,'color','w');
set(gca,'FontSize',24);
set(gca,'FontName','Arial');

%{
for ii = 1:size(VolumeSlice,1);
VolumeSlice(ii,1).FaceColor = 'interp';
VolumeSlice(ii,1).EdgeColor = 'none';
end
%}

VolumeSlice(1,1).FaceColor = 'none';
VolumeSlice(1,1).EdgeColor = 'none';
VolumeSlice(2,1).FaceColor = 'interp';
VolumeSlice(2,1).EdgeColor = 'none';
VolumeSlice(3,1).FaceColor = 'none';
VolumeSlice(3,1).EdgeColor = 'none';
hold on
theta = meshgrid(linspace(0, 2*pi, 50), linspace(0, 2*pi, 50)) ;

newX = a .* cos(theta);
newY = a .* sin(theta);

newZ = meshgrid(linspace(-L/2, L/2, 50), linspace(-L/2, L/2, 50))';

surf(newX,newY,newZ,'FaceColor','k','FaceAlpha',0.2,'LineStyle','none')
view([0 0])

%{
figure
colormap(hot)

%// Define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:size(PhiE_cavity,2), 1:size(PhiE_cavity,1));

%// Define a finer grid of points
[X2,Y2] = meshgrid(1:0.01:size(PhiE_cavity,2), 1:0.01:size(PhiE_cavity,1));

%// Interpolate the data and show the output
outData = interp2(X, Y, PhiE_cavity, X2, Y2, 'cubic');
imagesc(outData);

%// Cosmetic changes for the axes

set(gca, 'XTick', linspace(1,size(X2,2),11)); 
set(gca, 'YTick', linspace(1,size(X2,1),11));
set(gca, 'XTickLabel', linspace(-0.5,0.5,11));
set(gca, 'YTickLabel', linspace(-0.5,0.5,11));

%// Add colour bar
colorbar;

%imagesc(PhiE_cavity)
%{
sx = 0;
sy =0;
sz = z_z(length(z_z)/2);
VolumeSlice = slice(X,Y,Z,abs(E_3D),sx,sy,sz);
%}

%direction = [0 1 0];
%rotate(VolumeSlice,direction,90)

colorbar('vertical')
%{
VolumeSlice(1,1).FaceColor = 'none';
VolumeSlice(1,1).EdgeColor = 'none';
VolumeSlice(2,1).FaceColor = 'none';
VolumeSlice(2,1).EdgeColor = 'none';
VolumeSlice(3,1).FaceColor = 'interp';
VolumeSlice(3,1).EdgeColor = 'none';
%}
%{
hold on

theta = meshgrid(linspace(0, 2*pi, 50), linspace(0, 2*pi, 50)) ;

newX = a .* cos(theta);
newY = a .* sin(theta);

newZ = meshgrid(linspace(-L/2, L/2, 50), linspace(-L/2, L/2, 50))';

surf(newX,newY,newZ,'FaceColor','k','FaceAlpha',0.2,'LineStyle','--','LineWidth',2)
%}
axis equal

 ylabel('y/\mum')
 xlabel('x/\mum')
set(gcf,'color','w');
set(gca,'FontSize',24);
set(gca,'FontName','Arial');
view([0 90])
%}


