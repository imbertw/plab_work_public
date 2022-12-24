%neff_ring_sweep.m: Finds the effective index of a range of widths for a specific radius ring.

%parameters: clambda, width range, rOut, xs, nEffGuess, savefilename

%%
tic
%%Parameters

%Wavelength

clambda = 1.55;

%Layer stack for GF 45RFSOI process

xs = 'gf45spclo_xs';  

%Effective index guesses

nEffGuess= 2.6274;

%Sweep parameters

widths = 3.0:-0.05:0.2;
rOut = 40.;

%Save file name

savefile = 'nEffGuess_silicon_CBand.mat';

%%
c = 299792458;           % m/s
%GF45SPCLO
SZ = [];
for i=1:length(widths)
    
    for j = 1:length(rOut)
    
        % Find the effective index
        
        [neff(i,j) newSZ] = effective_index_ring(widths(i), rOut(j), clambda, xs, nEffGuess);

            if isempty(SZ)==0
                oldSZ = SZ;
            else
                oldSZ = [];
            end

        SZ = [oldSZ newSZ];
            if i<length(widths)
                nEffGuess = neff(i,j);
            else
            end
        save(savefile,'-v7.3')
    end
end
toc 
t = toc/60