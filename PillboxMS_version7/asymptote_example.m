clc;
clear all;
h=@(x) (x>=3 & x<=7).*1+(3<x & x>7).*3;
x = linspace(0,10,100);
figure(1)
plot(x,h(x));