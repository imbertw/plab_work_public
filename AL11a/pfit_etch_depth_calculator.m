%pfit_etch_depth_calculator: fits etch depth vs. resonance shift curve for rOut=9.97um and width =
%700nm in 45RF to a polynomial. Takes input of resonance shift required in microns,
%and outputs etch depth needed in microns.

function [etchdepth] = pfit_etch_depth_calculator(resonance_shift)
p = extract_polynomial();
etchdepth = polyval(p, resonance_shift);
end

function [p] = extract_polynomial()
load('filter_resonances_vs_etch_depth_r9970nm_w700nm.mat')
no_etch_resonance_wavelength = resonance_wavelength(1);
resonance_wavelength_shift = resonance_wavelength - no_etch_resonance_wavelength;
p = polyfit(resonance_wavelength_shift,etch_depth,2);
end
