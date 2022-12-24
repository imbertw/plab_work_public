%PillboxMV
function [mdvolume]=PillboxMV(a,L,n1,n2,n3,n4,h1,h2,q1,q2,beta1,beta2)
B = besselj(1,h1*a)/besselk(1,q1*a);
W = cos(beta1*L/2)/exp(1i*beta2*L/2);
intepsE_core_cavity = a*n1^2*pi*(a*h1*besselj(0,a*h1)^2-2*besselj(0,a*h1)*besselj(1,a*h1)+a*h1*besselj(1,a*h1)^2)*(beta1*L+sin(beta1*L))/(2*beta1*h1);
intepsE_clad_cavity = B^2*n2^2*pi^(1.5)*meijerG([],[1.5],[0 0 2],[],a^2*q1^2)*(beta1*L+sin(beta1*L))/(4*beta1*q1^2);
intepsE_core_mirror = 1i*a*exp(1i*beta2*L)*n3^2*pi*W^2*(a*h2*besselj(0,a*h2)^2-2*besselj(0,a*h2)*besselj(1,a*h2)+a*h2*besselj(1,a*h2)^2)/(beta2*h2);
intepsE_clad_mirror = 1i*B^2*exp(1i*beta2*L)*n4^2*pi^(1.5)*W^2*meijerG([],[1.5],[0 0 2],[],a^2*q2^2)/(2*beta2*q2^2);
maxEsq = n1^2*besselj(1,1.84118)^2;
mdvolume = (intepsE_core_cavity+intepsE_clad_cavity+intepsE_core_mirror+intepsE_clad_mirror)/maxEsq;
end
