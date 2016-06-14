clear all
close all
clc

Kb = 8.6173324*10^-5; %eV K^-1
Kb2 = 1.38064852*10^-23 
v1 = 10^13; %s^-1
beta = 1;  %K/s
fireamu = 6.64215616*10^-27 %kg

Tmax = [724.8, 698.4, 694.0, 709.1, 739.1, 753.4]';


redhead = @(T) Kb*T*(log(v1*T/beta)-3.64);


for n=1:length(Tmax)
    redhead(Tmax(n))
end


yes = @(p)  p*1/sqrt(2*pi*fireamu*Kb2*300);

10^15/0.0001/yes(1*10^-9*100000)
10^15/0.0001/yes(1*10^-12*100000)/60

10^15/(2.635*10^22*1*10^-9/(sqrt(2*300)))/60

kD2 = 2.330 * 10^-10 %1300 