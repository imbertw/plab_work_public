%Cylindrical_Parameter.m
function [h, q, beta] = CylindricalParameter(a, n1, n2, lambda)
error = 1e-2; %If we get the error message Error using fzero (line 257)Function values at interval endpoints must be finite and real, good chance we need to increase the error
k0 = 2*pi/lambda; %microns
betasq_TE01limit = n1^2*k0^2-(2.405/a)^2; %microns
betasq_TE02limit = n1^2*k0^2-(5.520/a)^2; %microns
betasq_min = n2^2*k0^2; %um^-2
betasq_max = n1^2*k0^2; %um^-2
betasq0 = [max([betasq_TE02limit betasq_min])+error min([betasq_TE01limit betasq_max])-error]; %need to add or subtract some error (e.g. 1e-2) just so we don't get infinities in the transcendental equation (neither h nor q can ever be zero)
betasq = fzero(@(beta2) TEcyl(beta2,n1,n2,k0,a), betasq0);
beta = sqrt(betasq); %um^-1
h = sqrt(n1^2*k0^2-beta^2); %um^-1
q = sqrt(beta^2-n2^2*k0^2); %um^-1
end

function [te0,h,q] = TEcyl(beta2,n1,n2,k0,a)
h=sqrt((n1*k0)^2-beta2);
q=sqrt(beta2-(n2*k0)^2);
%the objective is to find the zeroes of te0 and tm0
te0 = real(besselj(1,h.*a)./(h.*a.*besselj(0,h.*a))+besselk(1,q.*a)./(q.*a.*besselk(0,q.*a)));
end
