% m2cylphiR1.m - Discrete cylindrical modesolver for the E_phi field
% component.
% Imbert Wang, 3.20.2022
%% 0. Globals
%-------------------------------------------------------------------------%
% Global constants to be used in the code:
% u0 represents the magnetic permeability of free space, c is the speed of
% light, and Zo is the free-space impedance.
%-------------------------------------------------------------------------%

% c = 299792458;  % m/s
% Zo = u0*c;      % Impedance of free space
% e0 = 1/(u0*c^2);
micro = 1e-6;   %units for microns in meters
u0 = 4e-7*pi;   % H/m
%Set variables list
Dszn = round(logspace(1,3,5));
%% 1. Coordinate System and Pillbox radius
%-------------------------------------------------------------------------%
% First, we define our coordinate system.
%-------------------------------------------------------------------------%
rho_max_init = 1*micro;          %Maximum radius, units in metres
x_max_init   = 0.5*rho_max_init^2;   %Maximum x
a       = 0.2514*micro;          %Radius of core 0.2514um
%% 2. Structural and physical variables
%-------------------------------------------------------------------------%
%
%-------------------------------------------------------------------------%
lambda = 1.55*micro;             %Wavelength in meters
n1 = sqrt(12.085);               %Core index
n2 = sqrt(2.085);                %Cladding index
n3 = n2;                         %Core index in mirror layer
n4 = PillboxIndex(n1,n2);        %Cladding index in mirror layer


%% 3. Analytical pillbox modes
%Find transverse E modes
%n1 n2: cavity
[beta1,neff1,h1,q1,A1,B1,PhiE_cavity,Hr_cavity,Hz_cavity,r_r]=PillboxTV(a, rho_max_init, lambda, n1,n2); %run PillboxTV function to find transverse fields for n1/n2 layer
%n2 n3: mirror
[beta2,neff2,h2,q2,xx,xx,PhiE_mirror,Hr_mirror,Hz_mirror,r_r]=PillboxTV(a, rho_max_init, lambda, n3,n4); %run PillboxTV function to find transverse fields for n2/n3 layer
clear xx
%Generate longitudinal mode
[L,W,ZE1,ZE2,ZE3,z_z] = PillBoxLT(0,beta1,beta2); %run PillBoxLT function to get longitudinal mode
%Values to be passed to discretized modesolver
kz = beta1;
k0 = 2*pi/lambda;
n1eff = sqrt(n1^2 - (kz/k0)^2);
n2eff = sqrt(n2^2 - (kz/k0)^2);
%Loop starts here:
DX = zeros(1,length(Dszn));
LAMBDA_VEC = zeros(1,length(Dszn));
for ii=1:length(Dszn)
    N_core  = Dszn(ii);                   %Number of points for core
    dx      = 0.5*a^2/N_core;            %Discretization for the core
    DX(ii)  = dx;                        %Discretization saved variable
    N       = round(x_max_init/dx);  %Number of points
    x_max   = N*dx;                  %Maximum x-coordinate
    xvec    = 0:dx:x_max;           %x-coordinate, where x = rho^2
    rhovec  = (2*xvec).^(1/2);           %Radial coordinate
    rho_max = rhovec(end);           %Maximum radius for simulation
    
    
    %% 4. Setup for discrete modesolver
    nvec = zeros(N+1,1);
    nvec(rhovec <= a) =  n1eff;
    nvec(rhovec > a) = n2eff;
    %% 5. Wave Equation Operator
    %-------------------------------------------------------------------------%
    %
    %-------------------------------------------------------------------------%
    Hsqrtx = diag(xvec.^(1/2)); %Diagonal matrix of x vector
    
    % Create tridiagonal matrix
    Hdx2 = sparse(N+1,N+1);
    d = ones(N+1,1);
    Hdx2 = spdiags(d, -1, Hdx2);
    Hdx2 = spdiags(-2*d, 0, Hdx2);
    Hdx2 = spdiags(d, +1, Hdx2);
    Hdx2 = Hdx2/dx^2;
    % Diagonal k-vector
    %Ksq  = diag((k0*nvec).^2 - kz^2);
    H_inverseN = spdiags(1./nvec,0,N+1,N+1);
    Htot = -2*H_inverseN*Hsqrtx*Hdx2*Hsqrtx*H_inverseN;
    
    %% 6. Compute the fields and eigenvalues
    [V,D] = eigs(Htot,6,k0^2);
    %e_phi = V(:,1)./nvec;
    %e_phi = e_phi.'*sqrt(u0);
    lambda_vec = 2*pi./sqrt(diag(D));
    LAMBDA_VEC(ii) = lambda_vec(1);
    relerror(ii) = abs((LAMBDA_VEC(ii)-lambda)/lambda);
    clear Hsqrtx Hdx2 H_inverseN Htot nvec V D
end

%% 7. Plot wavelength eigenvalue vs. discretization

% figure;
% plot(DX*1e6,LAMBDA_VEC*1e9,'LineWidth',2)
% hold on
% plot(xlim,[lambda lambda]*1e9,'--','LineWidth',2)
% legend('discretized','analytical')
% xlabel('discretization /\mum')
% ylabel('wavelength eigenvalue /nm')
% title('Wavelength convergence vs. discretization of solver')

figure;
loglog(DX*1e9,relerror,'--o','LineWidth',2)
hold on
loglog(DX*1e9, DX*1e17*1e-4)
hold on
loglog(DX*1e9, (DX*1e17).^2*1e-4)
legend('relative error','O(x)', 'O(x^2')
%hold on
%plot(xlim,[lambda lambda]*1e9,'--','LineWidth',2)
%legend('discretized','analytical')
xlabel('discretization /(0.5*nm^2)')
ylabel('Wavelength error = |(\lambda_{discretized} - \lambda_{analytical})/\lambda_{analytical}|')
title('Wavelength convergence vs. discretization of solver')

% %% 7. Plot radial field
% rhovec_plot = unique([-fliplr(rhovec) rhovec]);
% E2_phi_plot = abs(e_phi).^2;
% max_E_phi_plot = max(E2_phi_plot);
% E2_phi_plot = [fliplr(E2_phi_plot(rhovec>0)) E2_phi_plot];
% figure;
% E2_analytical = PhiE_cavity(:,ceil(end/2)).^2;
% maxE2_analytical = max(E2_analytical);
% E2_analytical_plot = E2_analytical*max_E_phi_plot/maxE2_analytical;
% plot(r_r/micro,E2_analytical_plot,'LineWidth',2)
% hold on
% ylabel('|E_{\phi}|^2')
% xlabel('\rho/\mum')
% set(gcf,'color','w');
% set(gca,'FontSize',12);
% set(gca,'FontName','Arial');
% plot(rhovec_plot/micro,E2_phi_plot,'--','LineWidth',2)
% title(['Profile of fundamental TE mode for pillbox, a = ', num2str(a)])
% legend('analytical','discretized')
%eigenvalue = D(1,1);
% omega = sqrt(D(1,1)/(u0*e0));
% wavelength_calc = 2*pi*c/omega;