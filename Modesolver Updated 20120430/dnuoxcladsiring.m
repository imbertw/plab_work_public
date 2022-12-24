function [dnu1550 field t gammareturn] = dnuoxcladsiring(width, height, rOut, clambda, dxy, OPTS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Start clock
tic

%Opts defaults
if nargin < 4
    OPTS = struct;
end 
%Default to no plotting
if isfield(OPTS, 'plot') == 0
    OPTS.plot = 0;
end

if isfield(OPTS, 'plotf') == 0
    OPTS.plotf = 0;
end

%Default if no gamma guess is made
if isfield(OPTS, 'mu_guess') == 0
   OPTS.mu_guess = 3*2*pi/clambda*rOut;
end

%Default to adjacent FSRs
if isfield(OPTS, 'numFSRs') == 0
   OPTS.numFSRs = 1;
end

%Add modesolver path on Z: drive (might be faster to call locally)
if ispc
    addpath('Z:\Users\Cale\Code\Modesolver Updated 20120430');
elseif isunix
    addpath('/pbox/Users/Cale/Code/Modesolver Updated 20120430/');
end
%Constants
scrsz = get(0,'ScreenSize');


%Set up dielectric size matrix
nlyrs = ones(3, 3); 
side = 1.5;
top =1.5; 
bot = 1.5;



OPTS.NMODES_CALC =1;                    %Number of modes to calculate

%PML Parameters
OPTS.PMLwidth = [0 1 0 0];         
OPTS.PMLsigma = [1 1];
%Bent waveguide
OPTS.radius = rOut-width-side; 

dlyrsx = [side, width, side+side];     %Thickness of layers
dlyrsy = [top, height, bot];           %Height of layers


%Estimate FSR
tempn= [1.44 1.44 1.44; 1.44 index_Silicon(clambda) 1.44; 1.44 1.44 1.44];
tempk01 = 2*pi/(clambda+0.001);
tempk02 = 2*pi/(clambda-0.001);
[N,F,V] = sisolver3d2(tempn, dlyrsx, dlyrsy, dxy, tempk01, OPTS);
gamma1 = real(F.beta(1));
[N,F,V] = sisolver3d2(tempn, dlyrsx, dlyrsy, dxy, tempk02, OPTS);
gamma2 = real(F.beta(1));
fsr = round((0.002/(gamma2-gamma1))/0.001)*0.001;


blambda = clambda-OPTS.numFSRs*fsr;  %Lowest wavelength
ulambda = clambda+OPTS.numFSRs*fsr;  %Highest wavelength
dl = OPTS.numFSRs*fsr/10.;           %Discretization of wavelengths
%lambdas = [(blambda-dl),clambda,(ulambda+dl)];
%lambdas = [(blambda-2*dl):dl:(blambda+2*dl),clambda,(ulambda-2*dl):dl:(ulambda+2*dl)];

%Center wavelength plus 3 near each signal and idler
lambdas = [(blambda-dl):dl:(blambda+dl),clambda,(ulambda-dl):dl:(ulambda+dl)];


%Initialize vectors
gamma = zeros(1,length(lambdas));       %Initialize gamma
Sil = zeros(1,length(lambdas));         %Initialize index vs wavelength Si
SiO2array = zeros(1,length(lambdas));   %Initialize index vs wavelength SiO2
field = cell(2,1); 



for j = 1:length(lambdas)
    %Sellmeier equation for SiO2
    SiO2array(j) =  sqrt(1 + 0.6961663*(lambdas(j))^2/((lambdas(j))^2-0.0684043^2)+0.4079426*(lambdas(j))^2/((lambdas(j))^2-0.1162414^2)+0.8974794*(lambdas(j))^2/((lambdas(j))^2-9.896161^2));
    %NSiO2 = 1.44;
    NSiO2 = SiO2array(j);                 %Grab current index for SiO2
    NSi = index_Silicon(lambdas(j), 0);

    %Set up index layers
    nlyrs(1,1:3)=NSiO2;
    nlyrs(2,2) = NSi; nlyrs(2,1) = NSiO2; nlyrs(2,3) = NSiO2;
    nlyrs(3,1:3)=NSiO2; 

    k0 = 2*pi/lambdas(j);         %Input wavelength in units of wavenumber
    %Run modesolver
    [N,F,V] = sisolver3d2(nlyrs, dlyrsx, dlyrsy, dxy, k0, OPTS);

    %Store field profile 
    field{j,1} = real(F.Ex(:,:,1))';

    %Plot first mode
    if OPTS.plotf == 2
        m = 1;   %This is unnecessary but useful to add in a for loop if want to look at higher order modes
        figure;
        subplot(2,2,1); imagesc(N.x, N.y,(real(F.Ex(:,:,m))')); xlabel('y [micron]'); ylabel('x [micron]');
        title('Real Ex'); colormap(redbluehilight); colorbar; caxis([-2.5e-3,2.5e-3]);
        subplot(2,2,2); imagesc(N.x, N.y,(real(F.Ey(:,:,m))'));xlabel('y [micron]'); ylabel('x [micron]');
        title('Real Ey'); colormap(redbluehilight); colorbar; caxis([-2.5e-3,2.5e-3]);
        subplot(2,2,3); imagesc(N.x, N.y,(imag(F.Ez(:,:,m))'));xlabel('y [micron]'); ylabel('x [micron]');
        title('Imaginary Ez'); colormap(redbluehilight); colorbar; caxis([-2.5e-3,2.5e-3]);
    end
    %}
    gamma(1,j) = F.beta(1,1);                     %azimuthal propagation constant gamma
    OPTS.mu_guess = real(F.beta(1,1));            %Store previous gamma to use as guess for next
end

%Convert wavelengths to frequencies
omega = 2*pi./(lambdas)*2.99792458*10^8*10^6; %[1/s]

%Just plot the first mode
if OPTS.plotf == 1
    figure; imagesc(field{1,1}); 
end



%Also plot omega vs gamma, gamma+1 and gamma-1
if OPTS.plot == 1
    figure('Position',[10 50 0.99*scrsz(3) 0.89*scrsz(4)],'Renderer','zbuffer'); p3 = scatter(real(gamma), real(omega));  
    set(p3,'LineWidth', 2.5);
    hold all
    p3 = scatter(real(gamma+OPTS.numFSRs), real(omega));  
    set(p3,'LineWidth', 2.5);
    p3 = scatter(real(gamma-OPTS.numFSRs), real(omega));  
    set(p3,'LineWidth', 2.5);
    set(gca,'FontSize', 14); set(gcf, 'Color', [1 1 1]); 
    ylabel('\omega (rad/s)','FontSize', 16); xlabel('Gamma', 'FontSize', 16); title(['Oxide cladded Si waveguide (rOut = ', num2str(rOut), '\mum  w = ', num2str(width*1e3), 'nm  h = ', num2str(height*1e3), 'nm)'], 'FontSize', 20);
    grid on
    hold off
end
%}


%resample the lines
dg = 0.01; 
%+/- 2 instead of 1
xgamma = (min(real(gamma(1,:)))):dg:max(real(gamma(1,:)));
pomega = spline(real(gamma), real(omega),xgamma);
somega = spline(real(gamma+OPTS.numFSRs), real(omega),xgamma);
iomega = spline(real(gamma-OPTS.numFSRs), real(omega),xgamma);


if OPTS.plot ==1 
    %Plot resampled points
    figure('Position',[10 50 0.99*scrsz(3) 0.89*scrsz(4)],'Renderer','zbuffer'); p4 = scatter(xgamma, pomega);  
    set(p4,'LineWidth', 2.5);
    hold all
    p4 = scatter(xgamma, somega);  
    set(p4,'LineWidth', 2.5);
    p4 = scatter(xgamma, iomega);  
    set(p4,'LineWidth', 2.5);
    set(gca,'FontSize', 14); set(gcf, 'Color', [1 1 1]); 
    ylabel('\omega (rad/s)','FontSize', 16); xlabel('Gamma', 'FontSize', 16); title(['Oxide cladded Si waveguide (rOut = ', num2str(rOut), '\mum  w = ', num2str(width*1e3), 'nm  h = ', num2str(height*1e3), 'nm)'], 'FontSize', 20);
    grid on
    hold off
end


%Plot FSR difference vs gamma
dnu = (2*pomega-somega-iomega)/(2*pi*1e9);     %[GHz] 
plambda = 2*pi./pomega*2.99792458*10^8*10^9;
dnu1550 = spline(plambda,dnu,clambda*10^3);   %FSR mismatch at 1550


ddnudnu = diff(dnu)./diff(pomega/(2*pi*1e9)); 
difflambda = (plambda(2:end)+plambda(1:end-1))/2;
ddnudnu1550 = spline(difflambda,ddnudnu,clambda*10^3); 

if OPTS.plot == 1
    figure; plot(plambda, dnu); title('\Delta\nu vs wavelength'); ylabel('GHz'); 
    figure; plot(difflambda, ddnudnu); title('d\Delta\nu / d\nu vs wavelength'); ylabel('Hz/Hz'); 
    figure('Position',[10 50 0.99*scrsz(3) 0.89*scrsz(4)],'Renderer','zbuffer'); 
    [AX, H1, H2] = plotyy(xgamma, dnu,xgamma,plambda);  
    set(H1,'LineWidth', 2.5)
    set(H2,'LineStyle','--','LineWidth', 2.5)
    set(get(AX(1),'Ylabel'),'String','\Delta\nu (GHz)','FontSize', 16) 
    set(get(AX(2),'Ylabel'),'String','Wavelength (nm)','FontSize', 16) 
    set(get(AX(1),'Xlabel'),'String','Gamma (Mode Order) of Pump','FontSize', 16) 
    set(gcf, 'Color', [1 1 1]); 
    set(AX(1),'FontSize', 14); set(AX(2),'FontSize', 14);
    ylimits1 = get(AX(1),'YLim'); yinc1 = (ylimits1(2)-ylimits1(1))/10;
    ylimits2 = get(AX(2),'YLim'); yinc2 = (ylimits2(2)-ylimits2(1))/10;
    xlimits = get(AX(1),'XLim');
    xinc = 0.2;
    set(AX(1),'XTick',[xlimits(1):xinc:xlimits(2)],...
    'YTick',[ylimits1(1):yinc1:ylimits1(2)])
    set(AX(2),'YTick',[ylimits2(1):yinc2:ylimits2(2)])
    set(gca,'FontSize', 14);
    title(['FSR Matching for FWM (rOut = ', num2str(rOut), '\mum  w = ', num2str(width*1e3), 'nm  h = ', num2str(height*1e3), 'nm)', ' with \Delta\nu = 2\nu_{\gamma} - \nu_{\gamma-1} - \nu_{\gamma+1}'], 'FontSize', 20);
    grid on
end
%}
t = toc;

gammareturn = gamma(1,round(length(gamma(1,:))/2)); 
end

