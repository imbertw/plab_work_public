% Silicon refractive index fits to expt'al data
% from D.F. Edwards, Silicon (Si) in Palik, Handbook of Optical Constants
% of Solids (1985), p.547-.
%
% Syntax:  nn = index_Silicon(lambda_um, fittype)
%
% Input:    lambda_um   [1-vector] wavelengths in microns; if blank, a plot is given.
%           fittype     [scalar] 0=default is Herzberger, 1=Sellmeier
% Output:   nn          [1-vector] refractive index
%
% Milos Popovic, Mar 31, 2006

function nn = index_Silicon(lam, fittype)
if (nargin < 1) lam = [1.12:0.01:2]; end
if (nargin < 2) fittype = 0; end

% Herzberger fit
%A = 3.41906; B = 1.23172E-1; C = 2.65456E-2;
%D = -2.66511E-8; E = 5.45852E-14; L = 1./(lam.^2 - 0.028);
%nnherz = A + B*L + C*L.^2 + D*lam.^2 + E*lam.^4;
fherz = inline(['3.41906 + 1.23172E-1 ./(lam.^2 - 0.028) + '...
        '2.65456E-2 ./(lam.^2 - 0.028).^2 + -2.66511E-8*lam.^2 + '...
        '5.45852E-14 * lam.^4'], 'lam');
nnherz = fherz(lam);

% Sellmeier fit
L1 = 1.1071;
fsell = inline('sqrt(11.6858 + 0.939816./lam.^2 + 8.10461E-3 * L1^2./(lam.^2 - L1^2))','lam','L1');
nnsell = fsell(lam, L1);

if (nargin < 1)
    % Exptal data from Edwards:
    lamdata = 1+[.2 .372 .4 .532 .6 .696 .8 1];
    ndata = 3+[.5193 .5007 .4876 .4784 .4710 .4644 .4578 .4490];
    figure; plot(lam, nnherz, lam, nnsell, lamdata, ndata, 'o-', 'MarkerSize', 2);
    xlabel('Wavelength (\mum)'); ylabel('Refractive index');
    title('Si refractive index (D.F. Edwards in Palik, Handbk of Opt. Const., 1985)');
    %figformat(3.5, 7, 'Times'); filestampplot(gcf, mfilename);
    legend('Herzberger fit','Sellmeier fit','Experiments');
    
    lam2 = [1.5:0.01:1.6].';
    fnfit = fit(lam2-1.55, (fherz(lam2)+fsell(lam2,L1))/2, 'poly2');
    figure; plot(lam2, [fherz(lam2) fsell(lam2,L1) fnfit(lam2-1.55)]);
    xlabel('Wavelength (\mum)'); ylabel('Refractive index');
    title('Si refractive index: telecom range quadratic fit');
    %figformat(3.5, 7, 'Times'); filestampplot(gcf, mfilename);
    legend('Herzberger fit (Edwards)','Sellmeier fit (Edwards)','Quadratic mean fit');
    xlim([1.5 1.6]);
    
    nn = fnfit;
else
    if(fittype == 0)
        nn = nnherz;
    else
        nn = nnsell;
    end
end
%nn
