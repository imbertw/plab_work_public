%optimetchdepth: uses an optimizer to determine the etch depth required to bring about a
%given resonance shift.

function [etch_depth] = optimetchdepth(resonance_shift)
fun = @(y) abs(resonance_shift - etchdepth_vs_shift_calculator(y));
x1 = 0.04;
x2 = 0.06;
options = optimset('TolX',1e-6);
etch_depth = fminbnd(fun,x1,x2,options);
end

function [resonance_wavelength_shift] = etchdepth_vs_shift_calculator(etch_depth)
box_thickness = 0.2-etch_depth;
[lambda_res xx gammas xx] = ring_resonances_box(0.7, 9.97, 1.55, 'gf45rfsoi_xs_box', 2.0, box_thickness);
resonance_wavelength = lambda_res(gammas==71);
no_etch_resonance_wavelength = 1.557977295667556;
resonance_wavelength_shift = resonance_wavelength - no_etch_resonance_wavelength;
end