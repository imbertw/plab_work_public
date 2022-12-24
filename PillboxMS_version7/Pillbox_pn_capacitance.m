%Pillbox_pn_capacitance.m
q = 1.6e-19; %C
e0=8.85e-12; %Fm^-1
er=12.085; 
Na=1e18; %cm^-3
Nd = Na; %cm^-3

Vt = 0.0295; %V
ni = 1.5e10; %cm^-3
A = 1.986e-9; %cm^2

phi = Vt*log(Nd*Na/ni^2);
Va = linspace(-0.5,phi,100); %V
x = sqrt(2*er*e0/q*(1/Nd+1/Na)*(phi-Va));
C=er*e0*A./x;

semilogy(Va,C,'LineWidth',2)
grid on
