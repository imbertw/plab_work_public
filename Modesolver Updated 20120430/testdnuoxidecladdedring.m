%test
clear;

OPTS.plotf = 1;
OPTS.plot = 1;
OPTS.numFSRs = 1;

radii = 4.9:.1:5.1;

for j = 1:length(radii)
[dnu1550(j) field t gamma(j)] = dnuoxcladsiring(0.425, 0.220, radii(j), 1.55, [0.01 0.01], OPTS);
end 


figure; plot(radii,dnu1550,radii,gamma); grid on; set(gcf, 'Color', [1 1 1]); 
xlabel('Radius (\mm)'); ylabel('\Delta\nu \gamma');