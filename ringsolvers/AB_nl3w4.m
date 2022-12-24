% calculate the integration over the nonlinear material region for single
% ring cavity
% input: F, a 4x1 struct matrix, stores the fields; 
%        region_nl: nonliner material region
%        nomodes: which modes in F to calcualte; 
%        mode: 'c' linear cavity mode, or 'R' curved cavity mode; 
% output: A, effective field overlap for kai3_1111; 
%        [B1122, B1212, B1221], effective field overlap for kai3_1122, kai3_1212, kai3_1221; 

% notes: 
% 06/16/2011: the code assumes the field profile comes from ModeSolver; for arbitrary field profile, modification is needed; 
% 06/27/2011: updated AB coefficient taking into account of the fact that field direction is angle dependent for ring mode. for [100] wafer; for the triple-ring cavity; 
%
% To do: 
% - for general FWM, allow 4 different fields (F's) to calculate A,B coefficient, i.e., re-use function AB_nl3w(); 

% function [A, B1122, B1212, B1221] = AB_nl3w2(F, region_nl, nmodes, mode, Lcavity)
function [A, B] = AB_nl3w4(F, N, n_nl, nmodes, mode, Lcavity)

if (nargin < 5) 
    Lcavity = 1;   % unit assumed: um; for ring cavities, the radius information is already stored in the input variable "F"; 
    if mode ~= 'R'
    error('Please give the cavity length! \n'); 
    end; % for ring cavity, the input parameter Lcavity is not needed, since ring radius is implicit in F;  
end

Ix = N.nstruct;
Iy = N.nstruct;
Iz = N.nstruct;

Ixtrue = logical(Ix > n_nl-n_nl*1e-3 & Ix < n_nl+n_nl*1e-3);
Iytrue = logical(Iy > n_nl-n_nl*1e-3 & Iy < n_nl+n_nl*1e-3);
Iztrue = logical(Iz > n_nl-n_nl*1e-3 & Iz < n_nl+n_nl*1e-3);

Tx = imresize(Ixtrue,[length(F(1).Rr) length(F(1).Zr)], 'nearest');
Ty = imresize(Iytrue,[length(F(1).Rz) length(F(1).Zz)], 'nearest');
Tz = imresize(Iztrue,[length(F(1).Rz) length(F(1).Zr)], 'nearest');

% assignin('base', 'Ix',Ix)
% assignin('base','Iy', Iy)
% assignin('base', 'Iz', Iz)
% assignin('base', 'Tx', Tx)
% assignin('base', 'Ty', Ty)
% assignin('base', 'Tz', Tz)

for ii=1:4    % get the field in the nonlinear material region; ii denotes modes with different frequencies; 
    exfield(:,:,ii) = F(ii).Ex(:,:,nmodes(ii)).*Tx;
%     hyfield(:,:,ii) = F(ii).Hy(iexstart:iexend,jexstart:jexend,nmodes(ii));
    eyfield(:,:,ii) = F(ii).Ey(:,:,nmodes(ii)).*Ty;
%     hxfield(:,:,ii) = F(ii).Ey(ieystart:ieyend,jeystart:jeyend,nmodes(ii));
    ezfield(:,:,ii) = F(ii).Ez(:,:,nmodes(ii)).*Tz;
%     hzfield(:,:,ii) = F(ii).Hz(ihzstart:ihzend,jhzstart:jhzend,nmodes(ii));
    
end

% manipulate the field matrix above to have the same dimension so that the
% resulting fields are defined at the same spot. 
exfield = (exfield(1:end-1, :, :)+exfield(2:end, :, :) )/2; 
% hyfield = (hyfield(:, 1:end-1, :)+hyfield(:, 2:end, :) )/2; 
eyfield = (eyfield(:,1:end-1, :)+eyfield(:,2:end, :) )/2; 
% hxfield = (hxfield(1:end-1, :, :)+hxfield(2:end, :, :) )/2; 
ezfield = (ezfield(:, :, :)+ezfield(:, :, :) )/2; 

% for the FWM in triple-ring cavity, the uncoupld ring cavity mode will be
% used here(F(ii), ii=1,2,3,4 are the same)
exfield = exfield(:,:,1);
eyfield = eyfield(:,:,1);
ezfield = ezfield(:,:,1);
% assignin('base', 'exfield', exfield)
% assignin('base', 'eyfield', eyfield)
% assignin('base', 'ezfield', ezfield)

if mode=='C'
    fprintf('solving nonlinear coefficient in straight waveguide: \n'); 
    Rmtx = Lcavity.* ones(size(exfield,1),size(exfield,2)); % radius matrix (multiplied with cavity length) for use in integration with cartesian coordinate; 
else if mode =='R'
        fprintf('solving nonlinear coefficient in ring cavity: \n'); 
        RzNL = F(1).Rz;
        Rmtx = 2*pi.* RzNL(:)*ones(1, size(exfield,2));   % radius matrix (multiplied with cavity length) for use in integration with cylindrical coordinate; 
    else
        error('please choose eith ''b''-- straight waveguide mode, or ''w''-- cavity mode!\n'); 
    end
end

%assignin('base', 'Rmtx', Rmtx)

%% the following expression are specific for the FWM in single ring cavity. For the general expressions of AB
% coefficient (nonlinear coupling between uncoupoled cavity mode and anisotropic nonliear chi^3), please refer to function AB_nl3w.m

A = F(1).dx* F(1).dy* sum(sum((abs(eyfield).^4+ 3/4.*abs(exfield).^4+ 3/4.*abs(ezfield).^4+ abs(exfield).^2.*abs(ezfield).^2+ 1/2.* real(exfield.^2.*conj(ezfield.^2))) .*Rmtx)); 
B = F(1).dx* F(1).dy* sum(sum((3/4.*abs(exfield).^4+ 3/4.*abs(ezfield).^4+ abs(exfield).^2.*abs(ezfield).^2+ 4.*abs(exfield).^2.*abs(eyfield).^2+ 4.*abs(eyfield).^2.*abs(ezfield).^2+ 1/2.*real(conj(exfield).^2.*ezfield.^2)+ 2.*real(conj(eyfield).^2.*exfield.^2)+ 2.*real(conj(eyfield).^2.*ezfield.^2)) .*Rmtx)); 
% B = 1/4* F(1).dx* F(1).dy* sum(sum((3/4.*exfield.^4+ 3/4.*ezfield.^4- 1/2.*exfield.^2.*ezfield.^2+ 6.*exfield.^2.*eyfield.^2- 2.*eyfield.^2.*ezfield.^2) .*Rmtx)); 

% A = F(1).dx* F(1).dy* (sum(sum(conj(exfield(:,:,1)) .*exfield(:,:,2) .*exfield(:,:,3) .*conj(exfield(:,:,4) ) .*Rmtx)) ...
%   + sum(sum(conj(eyfield(:,:,1)) .*eyfield(:,:,2) .*eyfield(:,:,3) .*conj(eyfield(:,:,4) ) .*Rmtx)) ...
%   + sum(sum(conj(ezfield(:,:,1)) .*ezfield(:,:,2) .*ezfield(:,:,3) .*conj(ezfield(:,:,4) ) .*Rmtx))); 

% B1122 = F(1).dx* F(1).dy* (sum(sum(conj(exfield(:,:,1)) .*exfield(:,:,2) .*eyfield(:,:,3) .*conj(eyfield(:,:,4) ) .*Rmtx)) ...    %xxyy
%       + sum(sum(conj(exfield(:,:,1)) .*exfield(:,:,2) .*ezfield(:,:,3) .*conj(ezfield(:,:,4) ) .*Rmtx)) ...    %xxzz
%       + sum(sum(conj(eyfield(:,:,1)) .*eyfield(:,:,2) .*exfield(:,:,3) .*conj(exfield(:,:,4) ) .*Rmtx)) ...    %yyxx
%       + sum(sum(conj(eyfield(:,:,1)) .*eyfield(:,:,2) .*ezfield(:,:,3) .*conj(ezfield(:,:,4) ) .*Rmtx)) ...    %yyzz
%       + sum(sum(conj(ezfield(:,:,1)) .*ezfield(:,:,2) .*exfield(:,:,3) .*conj(exfield(:,:,4) ) .*Rmtx)) ...    %zzxx
%       + sum(sum(conj(ezfield(:,:,1)) .*ezfield(:,:,2) .*eyfield(:,:,3) .*conj(eyfield(:,:,4) ) .*Rmtx))) ;      %zzyy
% 
% B1212 = F(1).dx* F(1).dy* (sum(sum(conj(exfield(:,:,1)) .*eyfield(:,:,2) .*exfield(:,:,3) .*conj(eyfield(:,:,4) ) .*Rmtx)) ...    %xyxy
%       + sum(sum(conj(exfield(:,:,1)) .*ezfield(:,:,2) .*exfield(:,:,3) .*conj(ezfield(:,:,4) ) .*Rmtx)) ...    %xzxz
%       + sum(sum(conj(eyfield(:,:,1)) .*ezfield(:,:,2) .*eyfield(:,:,3) .*conj(ezfield(:,:,4) ) .*Rmtx)) ...    %yzyz
%       + sum(sum(conj(eyfield(:,:,1)) .*exfield(:,:,2) .*eyfield(:,:,3) .*conj(exfield(:,:,4) ) .*Rmtx)) ...    %yxyx
%       + sum(sum(conj(ezfield(:,:,1)) .*exfield(:,:,2) .*ezfield(:,:,3) .*conj(exfield(:,:,4) ) .*Rmtx)) ...    %zxzx
%       + sum(sum(conj(ezfield(:,:,1)) .*eyfield(:,:,2) .*ezfield(:,:,3) .*conj(eyfield(:,:,4) ) .*Rmtx))) ;      %zyzy
% 
% B1221 = F(1).dx* F(1).dy* (sum(sum(conj(exfield(:,:,1)) .*eyfield(:,:,2) .*eyfield(:,:,3) .*conj(exfield(:,:,4) ) .*Rmtx)) ...    %xyyx
%       + sum(sum(conj(exfield(:,:,1)) .*ezfield(:,:,2) .*ezfield(:,:,3) .*conj(exfield(:,:,4) ) .*Rmtx)) ...    %xzzx
%       + sum(sum(conj(eyfield(:,:,1)) .*exfield(:,:,2) .*exfield(:,:,3) .*conj(eyfield(:,:,4) ) .*Rmtx)) ...    %yxxy
%       + sum(sum(conj(eyfield(:,:,1)) .*ezfield(:,:,2) .*ezfield(:,:,3) .*conj(eyfield(:,:,4) ) .*Rmtx)) ...    %yzzy
%       + sum(sum(conj(ezfield(:,:,1)) .*exfield(:,:,2) .*exfield(:,:,3) .*conj(ezfield(:,:,4) ) .*Rmtx)) ...    %zxxz
%       + sum(sum(conj(ezfield(:,:,1)) .*eyfield(:,:,2) .*eyfield(:,:,3) .*conj(ezfield(:,:,4) ) .*Rmtx))) ;      %zyyz

 


