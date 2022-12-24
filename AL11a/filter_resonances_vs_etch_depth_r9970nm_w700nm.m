%filter_resonances_vs_etch_depth_r9970nm_w700nm.m: This code uses the modesolver to
%calculate the ring resonances of the rOut=9.97um Width = 700nm filter ring
%with varying etch depths.

%% Calculate

etch_depth = 0:0.005:0.2;

for ii=1:length(etch_depth)
    box_thickness = 0.2-etch_depth(ii);
    [lambda_res omegas gammas SZ] = ring_resonances_box(0.7, 9.97, 1.55, 'gf45rfsoi_xs_box', 2.0, box_thickness);
    resonance_wavelength(ii) = lambda_res(gammas==71);
    fprintf('For etch depth = %2.3f um, resonance wavelength = %2.3f um\n', etch_depth(ii), resonance_wavelength(ii)) 
    clear gammas lambda_res omegas SZ
end

save('filter_resonances_vs_etch_depth_r9970nm_w700nm.mat')

