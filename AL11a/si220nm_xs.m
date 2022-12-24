function[nlyrs,dlyrsx,dlyrsy,left_to_rout] = si220nm_xs(width,lambdas)
for ii=1:length(lambdas)
%Material parameters
nSi = index_Si_CMG(lambdas(ii));
nSiO2 = index_SiO2(lambdas(ii));
                    
%Layer refractive indices
nlyrs = ones(3, 3);
nlyrs(3,:) = nSiO2;
nlyrs(2,2) = nSi; nlyrs(2,[1 3]) = nSiO2; 
nlyrs(1,:) = nSiO2;


%Layer thickness and widths
wSi  = width;            %waveguide width
hSi  = 0.220;          %waveguide height
side = 3.0;
topOxide  = 3.0;
botOxide = 3.0;

dlyrsx = [side  wSi  side];         %microns wide            
dlyrsy = [botOxide hSi topOxide];           %microns tall
left_to_rout = side+wSi;
end