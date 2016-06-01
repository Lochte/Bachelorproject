% TPD-Temp calibration
clear all
close all
clc

TC = [350 380 410 440 470 500 530 560 590 620 650 680 710 740 770 800 830 860 890];

pyrometer = [321.3 322 322.8 323.2 323.5 324.5 326.6 331 339.7 353.7 371.2 394.7 417.7 443.9 468.7 492.7 516.2 539 563.6];

pyroKadd = pyrometer + 273.15;

plot(TC,pyroKadd)
hold on

kali = fit(TC(11:end)',pyroKadd(11:end)','poly1')
plot(kali)

korigeret = TC*kali.p1