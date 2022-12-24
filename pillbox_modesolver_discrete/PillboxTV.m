function [beta,neff,h,q,A,B,E_2D,Hr_2D,Hz_2D,r_r]=PillboxTV(radius,radius_limit,wavelength,n1,n2)
%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
% Inputs
a = radius; %microns
lambda = wavelength; %microns
k0 = 2*pi/lambda; %um^-1
%--------------------------------------------------------------------------
%Set up coordinate system
%--------------------------------------------------------------------------
N_points = 100;
factor = 2;
%r = linspace(0,factor*a,N_points);
r = linspace(0,radius_limit,N_points);
r1 = r(r<=a);
r2 = r(r>a);
r_r = [fliplr(-r2) fliplr(-r1) r1 r2];
RRR = r_r;
r_r = unique(r_r);
%for 3D
[X,Y] = meshgrid(r_r);
R = sqrt(X.^2 + Y.^2) + eps;
%--------------------------------------------------------------------------
%Obtain cylindrical parameters
%--------------------------------------------------------------------------
[h,q,beta] = CylindricalParameter(a,n1,n2,lambda); %beta in um^-1
neff = beta/k0; %unitless
%--------------------------------------------------------------------------
%Obtain Fields
%--------------------------------------------------------------------------
[E_phi,E_2D,Hr_2D,Hz_2D,A,B]=PillboxFields(h,q,beta,lambda,a,r1,r2,R);

%--------------------------------------------------------------------------
%Output to workspace
%--------------------------------------------------------------------------
%assignin('base','beta',beta)
%assignin('base','neff',neff)
%assignin('base','h',h)
%assignin('base','q',q)
%assignin('base','E_phi',E_phi)
%assignin('base','E_2D',E_2D)
%assignin('base','Hr_2D',Hr_2D)
%assignin('base','Hz_2D',Hz_2D)
%assignin('base','RRR',RRR)


%--------------------------------------------------------------------------
%Plot Fields
%--------------------------------------------------------------------------
%{
figure('NumberTitle', 'off', 'Name', 'E_phi field');
subplot(2,2,1)
plot(r_r,E_phi.^2,'b','LineWidth',3)
%}
%{
subplot(2,2,2)
surf(X,Y,(E_3D.^2),'EdgeColor','none')
view(90,90)
%}

end


