function []=ThreeDFieldPlotIsosurface(a,L,E_3D,X,Y,Z,z_z)

%--------------------------------------------------------------------------
%Plot Fields
%--------------------------------------------------------------------------

figure
%colormap(redbluehilight(1024))
s = isosurface(X,Y,Z,E_3D,0.8*max(max(max(E_3D))));
p = patch(s);
isonormals(X,Y,Z,E_3D,p)
view(3);
set(p,'FaceColor',[1.0 0 0]);
set(p,'EdgeColor','none');
camlight;
lighting gouraud;
hold on
s2 = isosurface(X,Y,Z,E_3D,0.8*min(min(min(E_3D))));
p2 = patch(s2);
isonormals(X,Y,Z,E_3D,p2)
view(3);
set(p2,'FaceColor',[0 0 1]);
set(p2,'EdgeColor','none');
%direction = [0 0 1];
%direction = [0 1 0];
%rotate(VolumeSlice,direction,90)
%colorbar('vertical');

%{
for ii = 1:size(VolumeSlice,1);
VolumeSlice(ii,1).FaceColor = 'interp';
VolumeSlice(ii,1).EdgeColor = 'none';
end
%}
% IsoSurface(1,1).FaceColor = 'none';
% IsoSurface(1,1).EdgeColor = 'none';
% IsoSurface(2,1).FaceColor = 'interp';
% IsoSurface(2,1).EdgeColor = 'none';
% IsoSurface(3,1).FaceColor = 'none';
% IsoSurface(3,1).EdgeColor = 'none';

hold on
theta = meshgrid(linspace(0, 2*pi, 50), linspace(0, 2*pi, 50)) ;

newX = a .* cos(theta);
newY = a .* sin(theta);

newZ = meshgrid(linspace(-L/2, L/2, 50), linspace(-L/2, L/2, 50))';

surf(newX,newY,newZ,'FaceColor','k','FaceAlpha',0.2,'LineStyle','none')