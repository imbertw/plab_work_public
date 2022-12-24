%dnusolver_check.m
%This code is a diagnostic tool for dnusolver. It generates most of the
%same outputs as dnusolver, but also generates plots and SZs to observe
%and compare.
%function [dFSR_GHz neff t Q groupIndex plambda xgamma SZ] = dnusolver_check(width, rOut, clambda, xs, nEffGuess)
function [dFSR_GHz neff t Q groupIndex FSR1_GHz FSR2_GHz Veff_um3 plambda xgamma SZ] = dnusolver_check(width, rOut, clambda, xs, nEffGuess)
tic
scrsz = get(0,'ScreenSize');
c = 299792458;           % m/s
um = 1e-6;               % micron conversion

%Set up vector of wavelengths to be checked

%% Wavelength range

%Estimate how low in wavelength you need to go
blambda = floor(1/(1/clambda+1/(3*pi*rOut))*(1/0.001))/(1/0.001);  %start wavelength rounded to 1nm
%Estimate how high in wavelength you need to go
ulambda = ceil(1/(1/clambda-1/(3*pi*rOut))*(1/0.001))/(1/0.001);   %end wavelength rounded to dlambda
lambdas = blambda:(ulambda-blambda)/4:ulambda;

% Calculate omega in Hz
omega = 2*pi./(lambdas*um)*c;

%% Mode number (azimuthal propagation constant)
gamma = zeros(1,length(lambdas));       %Initialize gamma

%% Set options for sisolver3d2
OPTS.NMODES_CALC = 4;                 %Number of modes you want to calculate
OPTS.eigmode = 'b';                   %beta mode
%PML parameters
OPTS.PMLwidth = [0 0.5 0 0]; % [left right bottom top]
OPTS.PMLsigma = [0.2 0.2];
%Effective index guess
k0 =2*pi/lambdas(1);
OPTS.mu_guess = nEffGuess*k0*rOut;     %guess gamma

OPTS.epsavg = 'simple'; % simple dielectric averaging at step changes
%%
SZ =[];
for ii = 1:length(lambdas)
    %Set up index layers
    [nlyrs, dlyrsx, dlyrsy, left_to_rout] = eval([xs '(width,lambdas(ii))']);
    
    rCenter = rOut-width/2.;
    
    %Calculate wavevector
    k0 = 2*pi/lambdas(ii); 
    
    %Initalize modesolver
    %dxy           = [0.02, 0.005];                   %Discretization
    dxy = [0.01, 0.0075];
    
    OPTS.radius   = (rOut-left_to_rout); % specify bent mode simulation and set left side of simulation domain
    %Run modesolver
    [N,F] = sisolver3d4(nlyrs, dlyrsx, dlyrsy, dxy, k0, OPTS);
    
    % Calculate the mode volume, but only for the central frequency
    if ii==3
        N1 = N;
        F1 = F;
        nlyrs1 = nlyrs;
        dlyrsx1 = dlyrsx;
        dlyrsy1 = dlyrsy;
    end
    
    % Generate SZ to inspect
    newSZ = [];
    newSZ.N = N;
    newSZ.F = F;

    
    if isempty(SZ)==0
        oldSZ = SZ;
    else
        oldSZ = [];
    end
    
    SZ = [oldSZ newSZ];
    
    %Find the azimuthal propagation constant gamma
    gamma(1,ii)     = F.beta(1,1);
    OPTS.mu_guess   = real(F.beta(1,1));
    %Convert from angular propagation constant to linear propagation
    %constant
    beta(ii) = gamma(1,ii)/rCenter;       % beta is the longitudinal propagation constant
    
    %Plot first mode
    m = 1;   %This is unnecessary but useful to add in a for loop if want to look at higher order modes
    figure;
    subplot(2,2,1); imagesc(N.x, N.y,(real(F.Ex(:,:,m))')); xlabel('y [micron]'); ylabel('x [micron]');
    title('Real Ex'); colormap(redbluehilight); colorbar; caxis([-2.5e-3,2.5e-3]);
    subplot(2,2,2); imagesc(N.x, N.y,(real(F.Ey(:,:,m))'));xlabel('y [micron]'); ylabel('x [micron]');
    title('Real Ey'); colormap(redbluehilight); colorbar; caxis([-2.5e-3,2.5e-3]);
    subplot(2,2,3); imagesc(N.x, N.y,(imag(F.Ez(:,:,m))'));xlabel('y [micron]'); ylabel('x [micron]');
    title('Imaginary Ez'); colormap(redbluehilight); colorbar; caxis([-2.5e-3,2.5e-3]);


end

%% Plotting
% Gamma vs omega
figure('Position',[10 50 0.99*scrsz(3) 0.89*scrsz(4)],'Renderer','zbuffer'); p = plot(real(omega), real(gamma));  
set(p,'LineWidth', 2.5);
set(gca,'FontSize', 14); set(gcf, 'Color', [1 1 1]); 
xlabel('\omega (rad/s)','FontSize', 16); ylabel('Gamma', 'FontSize', 16); title(['Oxide cladded Si waveguide (rOut = ', num2str(rOut), '\mum  w = ', num2str(width*1e3), 'nm  process: ', xs, ')'], 'FontSize', 20);
grid on

% Flip axes (omega vs gamma)
figure('Position',[10 50 0.99*scrsz(3) 0.89*scrsz(4)],'Renderer','zbuffer'); p2 = plot(real(gamma), real(omega));  
set(p2,'LineWidth', 2.5);
set(gca,'FontSize', 14); set(gcf, 'Color', [1 1 1]); 
ylabel('\omega (rad/s)','FontSize', 16); xlabel('Gamma', 'FontSize', 16); title(['Oxide cladded Si waveguide (rOut = ', num2str(rOut), '\mum  w = ', num2str(width*1e3), 'nm  process: ', xs, ')'], 'FontSize', 20);
grid on

figure('Position',[10 50 0.99*scrsz(3) 0.89*scrsz(4)],'Renderer','zbuffer'); p3 = scatter(real(gamma), real(omega));  
set(p3,'LineWidth', 2.5);
hold all
p3 = scatter(real(gamma+1), real(omega));  
set(p3,'LineWidth', 2.5);
p3 = scatter(real(gamma-1), real(omega));  
set(p3,'LineWidth', 2.5);
set(gca,'FontSize', 14); set(gcf, 'Color', [1 1 1]); 
ylabel('\omega (rad/s)','FontSize', 16); xlabel('Gamma', 'FontSize', 16); title(['Oxide cladded Si waveguide (rOut = ', num2str(rOut), '\mum  w = ', num2str(width*1e3), 'nm  process: ', xs, ')'], 'FontSize', 20);
grid on
hold off

%Plot the structure
figure; imagesc(N.x, N.y, (N.n).'); set(gca, 'YDir', 'normal'); colorbar;

%% Interpolate the angular frequencies for finer sampling
%resample the lines
dg = 0.01; 
xgamma = (min(real(gamma(1,:)))+1):dg:max(real(gamma(1,:))-1);
pomega = spline(real(gamma), real(omega),xgamma);                          %Pump frequencies
somega = spline(real(gamma+1), real(omega),xgamma);                        %Signal frequencies (where there is a pump, signal, and idler)
iomega = spline(real(gamma-1), real(omega),xgamma);                        %Idler frequencies (where there is a pump, signal, and idler)

%% Calculate FSR difference vs gamma
dnu      = (2*pomega-somega-iomega)/(2*pi*1e9);     %FSR mismatches in GHz
plambda  = 2*pi./pomega*c*1e9;                      %Pump wavelength in nm
dFSR_GHz = spline(plambda,dnu,clambda*1e3);         %FSR mismatch at 1550

%% Calculate the FSR
fsr1     = (pomega-somega)/(2*pi*1e9);       %FSRs in GHz
fsr2     = (iomega-pomega)/(2*pi*1e9); 
FSR1_GHz = spline(plambda,fsr1,clambda*1e3); %FSRs at 1550
FSR2_GHz = spline(plambda,fsr2,clambda*1e3);

% Plot the omega vs. gamma
figure('Position',[10 50 0.99*scrsz(3) 0.89*scrsz(4)],'Renderer','zbuffer'); p4 = scatter(xgamma, pomega);  
set(p4,'LineWidth', 2.5);
hold all
p4 = scatter(xgamma, somega);  
set(p4,'LineWidth', 2.5);
p4 = scatter(xgamma, iomega);  
set(p4,'LineWidth', 2.5);
set(gca,'FontSize', 14); set(gcf, 'Color', [1 1 1]); 
ylabel('\omega (rad/s)','FontSize', 16); xlabel('Gamma', 'FontSize', 16); title(['Oxide cladded Si waveguide (rOut = ', num2str(rOut), '\mum  w = ', num2str(width*1e3), 'nm  process:', xs, ')'], 'FontSize', 20);
grid on
hold off
%% Plot the FSR difference vs. gamma
figure; plot(plambda, dnu); title('\Delta\nu vs wavelength'); ylabel('GHz');  
figure('Position',[10 50 0.99*scrsz(3) 0.89*scrsz(4)],'Renderer','zbuffer'); 
[AX, H1, H2] = plotyy(xgamma, dFSR_GHz,xgamma,plambda);  
set(H1,'LineWidth', 2.5)
set(H2,'LineStyle','--','LineWidth', 2.5)
set(get(AX(1),'Ylabel'),'String','\Delta\nu (GHz)','FontSize', 16) 
set(get(AX(2),'Ylabel'),'String','Wavelength (nm)','FontSize', 16) 
set(get(AX(1),'Xlabel'),'String','Gamma (Mode Order) of Pump','FontSize', 16) 
set(gcf, 'Color', [1 1 1]); 
set(AX(1),'FontSize', 14); set(AX(2),'FontSize', 14);
ylimits1 = get(AX(1),'YLim'); yinc1 = (ylimits1(2))/10;
ylimits2 = get(AX(2),'YLim'); yinc2 = (ylimits2(2)-ylimits2(1))/10;
xlimits = get(AX(1),'XLim');
xinc = 1;
set(AX(1),'XTick',[xlimits(1):xinc:xlimits(2)],...
    'YTick',[ylimits1(1):yinc1:ylimits1(2)])
set(AX(2),'YTick',[ylimits2(1):yinc2:ylimits2(2)])
set(gca,'FontSize', 14);
title(['FSR Matching for FWM (rOut = ', num2str(rOut), '\mum  w = ', num2str(width*1e3), 'nm  process: ', xs, ')', ' with \Delta\nu = 2\nu_{\gamma} - \nu_{\gamma-1} - \nu_{\gamma+1}'], 'FontSize', 20);
grid on
%% Calculate the effective index
neff = real(beta(3))/(2*pi/clambda);   %Effective index at center wavelength
%% Calculate the mode volume of the ring
n_nl = index_Si_CMG(clambda); % index of silicon at this frequency
G = [F1, F1, F1, F1]; %just use the same mode, as signal pump and idler modes are all pretty similar
nmodes = [1 1 1 1];   % use fhe fundamental mode
%region_nl = [dlyrsx(1),dlyrsx(1)+dlyrsx(2),dlyrsy(1),dlyrsy(1)+dlyrsy(2)];
%[Veff_um3]   = nonlinear_mode_volume_old(N1, G, nmodes, region_nl, n_nl, clambda); %the nonlinear interaction volume in microns^3
[Veff_um3]   = nonlinear_mode_volume(N1, G, nmodes, n_nl, clambda); %the nonlinear interaction volume in microns^3
%% Calculate the group index of the ring
%resample omega
dOmega       = (max(omega)-min(omega))/100;
xOmega_plus  = omega(3) + dOmega;
xOmega_minus = omega(3) - dOmega;
%Find the gradient
beta_plus  = spline(omega,real(beta),xOmega_plus);
beta_minus = spline(omega,real(beta),xOmega_minus);

betas = [beta_minus real(beta(3)) beta_plus];

dBetadOmega      =  diff(betas*1e6)/dOmega;    %s/m

groupVelocity    = 1./dBetadOmega;             %m/s
avgGroupVelocity = mean(groupVelocity);
groupIndex       = c/avgGroupVelocity;            % group index n_g

%% Calculate the loss Q of the ring
alpha = imag(beta(3));
Q = k0*groupIndex/(2*alpha);          %Q is the quality factor

%% end time
t = toc;
end