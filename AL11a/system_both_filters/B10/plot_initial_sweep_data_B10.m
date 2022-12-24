load('system_1_thru_laser_sweep_laser_sweep_11_19_2022_16_0.mat');
plot(XData, YData);
hold on;
clear all;
load('system_1_ASE_bypass_in_signal_out_laser_sweep_11_19_2022_17_15.mat');
plot(XData, YData);
hold on;
load('system_1_ASE_bypass_in_idler_out_laser_sweep_11_21_2022_19_41.mat');
plot(XData, YData);
set(gcf,'Color','w')
title('B10 System Site 1 Pre-etch')
ylabel('Transmission /dB')
xlabel('Wavelength /nm')
