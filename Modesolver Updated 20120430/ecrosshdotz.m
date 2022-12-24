%---.----.----.----.----.----.----.----.----.----.----.----.----.----.----.----.
% ecrosshdotz - Vector field conjugated or non-conjugated cross-product integral
%                    I = int(ei x hj* dot z^, dA), or
%                    I = int(ei x hj dot z^, dA)
%             - enables compution of various coupled-mode theory normalizations
%               e.g. S_{ij} = (1/2) int(ei x hj* dot z^, dA), or
%                    S_{ij} = (1/4) int(ei x hj* + ej* x hi dot z^, dA)
%
% Syntax:   [Pz, SzTE, SzTM] = ecrosshdotz(Fi, Fj, [imode jmode], flagcnc, x, y);
%                                                  {---------OPTIONAL---------}
%
% Inputs:   Fi, Fj   - data structures which contain electric fields Ex(M-1,N-2), Ey(M-2,N-1) and
%                      magnetic fields Hx(M-2,N-1), Hy(M-1,N-2) plus spacings dx, dy for integration
%                      which must be the same for both fields! Also must have Pr, Pz for complex domains.
%           [modes]  - (OPTIONAL) 2-vector, if given specifies mode of Fi and Fj to use, e.g. [2 1]
%                      Default value is [1 1].
%           flagcnc  - (OPTIONAL) 0 for non-conjugated (default), 1 for conjugated overlap
%           x, y     - (OPTIONAL)
%
% Milos Popovic, Nov 26, 2002

% Updates:
% --------
% Nov 26, 2002 - ecrosshdotz first code written.
% Mar  8, 2004 - Added outputs for the TE (Ex) and TM (Ey) power densities
% Jul 16, 2004 - Updated for a conjugated and non-conjugated integral version.
% Sep 15, 2005 - Updated for non-uniform grids, passed in as x,y, vectors
%                typically from N index data structure of mode solver.

function [Pz,SzTE,SzTM] = ecrosshdotz(Fi, Fj, modes, flagcnc, x, y)
PEC = 0; PMCh = 1; PBC = 2; PMC = 3; PECh = 4;              % Define boundary condition identifiers [2006-11-08 MP; added]

MX = size(Fi.Ex,1)+1; MY = size(Fi.Ey,2)+1;                 % Extract grid size MX and MY = #pixels+1
if (nargin < 3) modes = [1 1]; end                          % DEFAULT values for mode indices
if (nargin < 4) flagcnc = 0; end                            % DEFAULT is non-conjugated overlap
if (length(modes) < 2) modes(2) = modes(1); end             % If second mode index not given assume same.
if (nargin < 5)
    if (isfield(Fi,'dx') && isfield(Fi,'dy'))  dATE = Fi.dx*Fi.dy; dATM = dATE;
    else  error('ECROSSHDOTZ.m - must pass in full coordinate vectors x,y; *or* Fi must have Fi.dx, Fi.dy for uniform grids.');  end
else                                                        % Non-uniform grid
    % Find extra pixels where needed
    iL = 0; iR = 0; iD = 0; iU = 0;
    if(isfield(Fi,'ver') && strcmp(Fi.ver,'R2'))              % Change size depending on boundary conditions
        iL = any(Fi.BC(1) == [PMC PMCh PBC]); %any(strcmp(F.BC(1),{'PMC','PMCh','PBC'}));
        iR = any(Fi.BC(2) == [PMC PMCh]);     %any(strcmp(F.BC(2),{'PMC','PMCh'}));
        iD = any(Fi.BC(3) == [PMC PMCh PBC]); %any(strcmp(F.BC(3),{'PMC','PMCh','PBC'}));
        iU = any(Fi.BC(4) == [PMC PMCh]);     %any(strcmp(F.BC(4),{'PMC','PMCh'}));
    end

%    dxtmp = (x(3:2:end)-x(1:2:end-2)); dytmp = (y(4:2:end-1)-y(2:2:end-3));
    dxtmp = x(3:2:end)-x(1:2:end-2); dytmp = y([4-2*iD:2:end-1, end*ones(1,iU)])-y([ones(1,iD), 2:2:end-3]);    % [2006-11-08] added support for new R2
    if(iD == PMCh) dytmp(1) = 2*dytmp(2); end;  if(iU == PMCh) dytmp(end) = 2*dytmp(end); end;
    dATE = dxtmp(:) * dytmp(:).';                           % TE grid pixel areas
%    dxtmp = (x(4:2:end-1)-x(2:2:end-3)); dytmp = (y(3:2:end)-y(1:2:end-2));
    dxtmp = x([4-2*iL:2:end-1, end*ones(1,iR)])-x([ones(1,iL), 2:2:end-3]); dytmp = y(3:2:end)-y(1:2:end-2);  % [2006-11-08] added support for new R2
    if(iL == PMCh) dxtmp(1) = 2*dxtmp(2); end;  if(iR == PMCh) dxtmp(end) = 2*dxtmp(end); end;
    dATM = dxtmp(:) * dytmp(:).';                           % TM grid pixel areas
end

mi = modes(1); mj = modes(2);

% Find (e* dot delta_epsilon dot e) integrand components
if(flagcnc == 1)
    SzTE = -Fi.Ex(:,:,mi) .* conj(Fj.Hy(:,:,mj));
    SzTM = Fi.Ey(:,:,mi) .* conj(Fj.Hx(:,:,mj));
else
    SzTE = -Fi.Ex(:,:,mi) .* (Fj.Hy(:,:,mj));
    SzTM = Fi.Ey(:,:,mi) .* (Fj.Hx(:,:,mj));
end

%Pz = sum([SzTE(:); SzTM(:)]) * Fi.dx*Fi.dy;                 % Perform integral
Pz = sum([SzTE(:).*dATE(:); SzTM(:).*dATM(:)]);             % Perform integral for non-uniform or uniform grid
% (minus sign because our fields are in fact in x and z coordinates, dot y)
