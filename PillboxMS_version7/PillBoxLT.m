function [L,W,ZE1,ZE2,ZE3,z_z] = PillboxLT(m,beta1,beta2)
L = 2/beta1*atan(-1i*beta2/beta1);
U = mod(m+1,2); V = mod(m,2);
if U==0
    pm = -1;
else
    pm = 1;
end
W = U*cos(beta1*L/2)/exp(1i*beta2*L/2)+V*sin(beta1*L/2)/exp(1i*beta2*L/2);
%z_z = linspace(-L,L,100);
z_z = linspace(-.5,.5,100);
z1 = z_z(z_z<=-L/2);
z2 = z_z(z_z>-L/2 & z_z<L/2);
z3 = z_z(z_z>=L/2);
ZE2 = U*cos(beta1*z2)+V*sin(beta1*z2);
ZE3 = W*exp(1i*beta2*z3);
ZE1 = pm*W*exp(-1i*beta2*z1);
end
