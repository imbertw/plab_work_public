% Compute group index from modesolver fields
% - assumes no material dispersion (for more general case, consult Memo 009, from Dec 2004)
%
% Syntax:   ng = groupindex(F, N, {modesvec})
%
% Input:    F - Mode data structure containing fields & discretization (Ex,Ey,Ez,Hx,Hy,Hz,dx,dy)
%           N - Data structure containing index distribution of waveguide
%           modesvec - [OPTIONAL] for which modes to find group index (default: all)
% Output:   ng - vector of group indices
%
% Milos Popovic, Apr 15, 2005

% Code updates:
% -------------
% Apr 15, 2005 - First version.
% Sep 11, 2005 - Updated for non-uniform grids.
% Aug  5, 2006 - Added line to ensure N.x, N.y are vertical vectors
% Nov  8, 2006 - Added ability to deal with solver engine R2 modes with
%                various boundary conditions (PEC,PMC,PECh,PMCh,PBC).

function ng = groupindex(F, N, modesvec)
PEC = 0; PMCh = 1; PBC = 2; PMC = 3; PECh = 4;                              % Define boundary condition identifiers [2006-11-08 MP; added]
c = 299792458; u0 = 4e-7*pi; e0 = 1/c^2/u0;
if (nargin < 3)  modesvec = [1:size(F.Ex,3)];  end;     % If not specified find group index of all modes

iL = 0; iR = 0; iD = 0; iU = 0;
if(isfield(F,'ver') && strcmp(F.ver,'R2'))              % Change size depending on boundary conditions
    iL = any(F.BC(1) == [PMC PMCh PBC]); %any(strcmp(F.BC(1),{'PMC','PMCh','PBC'}));
    iR = any(F.BC(2) == [PMC PMCh]);     %any(strcmp(F.BC(2),{'PMC','PMCh'}));
    iD = any(F.BC(3) == [PMC PMCh PBC]); %any(strcmp(F.BC(3),{'PMC','PMCh','PBC'}));
    iU = any(F.BC(4) == [PMC PMCh]);     %any(strcmp(F.BC(4),{'PMC','PMCh'}));
end

for k = 1:length(modesvec)
    mm = modesvec(k);

    if (isfield(N,'x') & isfield(N,'y'))
         N.x = N.x(:); N.y = N.y(:);        % [20060805] Added just to make sure orientation is correct
% %        I1 = sum(sum(conj(F.Ex(:,:,mm)) .* N.n(2:2:end-1,3:2:end-2).^2 .* F.Ex(:,:,mm))) + ...
% %             sum(sum(conj(F.Ey(:,:,mm)) .* N.n(3:2:end-2,2:2:end-1).^2 .* F.Ey(:,:,mm))) + ...
% %             sum(sum(conj(F.Ez(:,:,mm)) .* N.n(3:2:end-2,3:2:end-2).^2 .* F.Ez(:,:,mm)));
% %        I2 = sum(sum(conj(F.Hx(:,:,mm)) .* F.Hx(:,:,mm))) + ...
% %             sum(sum(conj(F.Hy(:,:,mm)) .* F.Hy(:,:,mm))) + ...
% %             sum(sum(conj(F.Hz(:,:,mm)) .* F.Hz(:,:,mm)));
% %        I3 = 2*real(ecrosshdotz(F, F, [mm mm], 1));
%         I1 = sum(sum(conj(F.Ex(:,:,mm)) .* N.n(2:2:end-1,3:2:end-2).^2 .* F.Ex(:,:,mm) .* ...
%                     ((N.x(3:2:end)-N.x(1:2:end-2))   * (N.y(4:2:end-1)-N.y(2:2:end-3)).') )) + ...      % Non-uniform grid pixel areas
%              sum(sum(conj(F.Ey(:,:,mm)) .* N.n(3:2:end-2,2:2:end-1).^2 .* F.Ey(:,:,mm) .* ...
%                     ((N.x(4:2:end-1)-N.x(2:2:end-3)) * (N.y(3:2:end)-N.y(1:2:end-2)).')   )) + ...
%              sum(sum(conj(F.Ez(:,:,mm)) .* N.n(3:2:end-2,3:2:end-2).^2 .* F.Ez(:,:,mm) .* ...
%                     ((N.x(4:2:end-1)-N.x(2:2:end-3)) * (N.y(4:2:end-1)-N.y(2:2:end-3)).') ));
%         I2 = sum(sum(conj(F.Hx(:,:,mm)) .* F.Hx(:,:,mm) .* ...
%                     ((N.x(4:2:end-1)-N.x(2:2:end-3)) * (N.y(3:2:end)-N.y(1:2:end-2)).')   )) + ...
%              sum(sum(conj(F.Hy(:,:,mm)) .* F.Hy(:,:,mm) .* ...
%                     ((N.x(3:2:end)-N.x(1:2:end-2))   * (N.y(4:2:end-1)-N.y(2:2:end-3)).') )) + ...
%              sum(sum(conj(F.Hz(:,:,mm)) .* F.Hz(:,:,mm) .* ...
%                     ((N.x(3:2:end)-N.x(1:2:end-2))   * (N.y(3:2:end)-N.y(1:2:end-2)).')   ));
%        cc = sqrt(2);
%        if(iL) F.Ey(1,:,m) = cc*F.Ey(1,:,m); F.Ez(1,:,m) = cc*F.Ez(1,:,m); F.Hx(1,:,m) = cc*F.Hx(1,:,m); end
%        if(iR) F.Ey(M-1+iL,:,m) = cc*F.Ey(M-1+iL,:,m); F.Ez(M-1+iL,:,m) = cc*F.Ez(M-1+iL,:,m); F.Hx(M-1+iL,:,m) = cc*F.Hx(M-1+iL,:,m); end
%        if(iD) F.Ex(:,1,m) = cc*F.Ex(:,1,m); F.Ez(:,1,m) = cc*F.Ez(:,1,m); F.Hy(:,1,m) = cc*F.Hy(:,1,m); end
%        if(iU) F.Ex(:,N-1+iD,m) = cc*F.Ex(:,N-1+iD,m); F.Ez(:,N-1+iD,m) = cc*F.Ez(:,N-1+iD,m); F.Hy(:,N-1+iD,m) = cc*F.Hy(:,N-1+iD,m); end

        dxr = (N.x(3:2:end)-N.x(1:2:end-2));
        dyr = (N.y([4-2*iD:2:end-1, end*ones(1,iU)])-N.y([ones(1,iD), 2:2:end-3]));
        dxz = (N.x([4-2*iL:2:end-1, end*ones(1,iR)])-N.x([ones(1,iL), 2:2:end-3]));
        dyz = (N.y(3:2:end)-N.y(1:2:end-2));
        
        I1 = sum(sum(conj(F.Ex(:,:,mm)) .* N.n(2:2:end-1,3-2*iD:2:end-2+2*iU).^2 .* F.Ex(:,:,mm) .* ...
                    (dxr   * dyr.') )) + ...      % Non-uniform grid pixel areas
             sum(sum(conj(F.Ey(:,:,mm)) .* N.n(3-2*iL:2:end-2+2*iR,2:2:end-1).^2 .* F.Ey(:,:,mm) .* ...
                    (dxz * dyz.')   )) + ...
             sum(sum(conj(F.Ez(:,:,mm)) .* N.n(3-2*iL:2:end-2+2*iR,3-2*iD:2:end-2+2*iU).^2 .* F.Ez(:,:,mm) .* ...
                    (dxz * dyr.') ));
        I2 = sum(sum(conj(F.Hx(:,:,mm)) .* F.Hx(:,:,mm) .* ...
                    (dxz * dyz.')   )) + ...
             sum(sum(conj(F.Hy(:,:,mm)) .* F.Hy(:,:,mm) .* ...
                    (dxr   * dyr.') )) + ...
             sum(sum(conj(F.Hz(:,:,mm)) .* F.Hz(:,:,mm) .* ...
                    (dxr   * dyz.')   ));
        I3 = 2*real(ecrosshdotz(F, F, [mm mm], 1, N.x, N.y));
    else
        error('GROUPINDEX.m - Error: must have N.x and N.y coordinates supplied for integration in call groupindex(F,N,mode).');
    end

%    ng(k,1) = c * (e0*I1 + u0*I2)*F.dx*F.dy / I3;
    ng(k,1) = c * (e0*I1 + u0*I2) / I3;
end
