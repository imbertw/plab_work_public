%plot_etch_depth_vs_resonance_shift_rOut_9_97um_width_700nm.m: plot the etch depth vs. resonance shift for rOut = 9.97um Width = 700nm
%45RF rings
load('filter_resonances_vs_etch_depth_r9970nm_w700nm.mat')
no_etch_resonance_wavelength = resonance_wavelength(1);
resonance_wavelength_shift = resonance_wavelength - no_etch_resonance_wavelength;


p = plot(etch_depth, resonance_wavelength_shift*1e3, '-o', 'LineWidth', 2)
grid on
grid minor
title('Filter ring resonances for different oxide etch depths, ROut=9970nm, Width=700nm')
xlabel('Etch Depth /\mum')
ylabel('Resonance Wavelength Shift /nm')
set(gcf,'Color','w')

