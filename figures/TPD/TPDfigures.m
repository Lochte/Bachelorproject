function TPD()
clear all
close all
clc

%% Kalibrering

kali = importdata('kalibrering.csv',';',4);

TC = kali.data(:,1);
pyroC = kali.data(:,2);
pyroK = kali.data(:,2)+273.15*ones(length(pyroC),1);

plot(TC,pyroK,'*')

kalifit = fit(TC(10:end),pyroK(10:end),'poly1');
hold on
plot(kalifit)
legend ('Datapoints','Linear fit','Location','NorthWest')

ylabel('Pyrometer temperature [K]')
xlabel('Thermocouple temperature [K]')
set(gca,'FontSize',14,'FontWeight','bold')

print(['kalibrering'],'-dpng','-r200')

%% Plots

figure
TPDplot('GrIr1740CD2',kalifit.p1,kalifit.p2)

figure
TPDplot('GrIr1740CD2longanneal',kalifit.p1,kalifit.p2)

figure
TPDplot('GrIr1740CD2_2504',kalifit.p1,kalifit.p2)

figure
TPDplottemp('GrIr1740CD2',kalifit.p1,kalifit.p2,'b',650)

figure
TPDplottemp('GrIr1740CD2_2504',kalifit.p1,kalifit.p2,'r',2)
hold on
TPDplottemp('GrIr1740CD2longanneal',kalifit.p1,kalifit.p2,'g',2)
legend({'Gr/Ir-1740C dose-D2' 'Gr/Ir-1740C dose-D2-after long anneal'},'Location','NorthWest')
print(['GrIr1740CD2longanneal' 'temp'],'-dpng','-r200')


figure
TPDplot('050516IrSorenD2dose1hour',kalifit.p1,kalifit.p2)

figure
TPDplot('050516IrSorenD2dose1hour2',kalifit.p1,kalifit.p2)

figure
TPDplot('050516IrSorenDdose1hour2',kalifit.p1,kalifit.p2)

figure
TPDplottemp('050516IrSorenD2dose1hour2',kalifit.p1,kalifit.p2,'b',2)
hold on 
TPDplottemp('050516IrSorenD2dose1hour',kalifit.p1,kalifit.p2,'g',2)
legend({'Gr/Ir-1740C dose-D2-#1' 'Gr/Ir-1740C dose-D2-#2'},'Location','NorthWest')
print(['050516IrSorenD2dose1hour' 'temp'],'-dpng','-r200')



figure
TPDplottemp('050516IrSorenDdose1hour2',kalifit.p1,kalifit.p2,'b',2)
hold on
TPDplottemp('050516IrSorenDdose1hour',kalifit.p1,kalifit.p2,'g',2)
legend({'Gr/Ir-1740C dose-D-#1' 'Gr/Ir-1740C dose-D-#2'},'Location','NorthWest')
print(['050516IrSorenDdose1hour' 'temp'],'-dpng','-r200')



figure
peakfit('050516IrSorenD2dose1hour2',kalifit.p1,kalifit.p2,'b',2,723,2)
figure
peakfit('050516IrSorenD2dose1hour',kalifit.p1,kalifit.p2,'g',2,700,1)
figure
peakfit('050516IrSorenDdose1hour2',kalifit.p1,kalifit.p2,'b',2,620,1)
figure
peakfit('050516IrSorenDdose1hour',kalifit.p1,kalifit.p2,'g',2,700,1)
figure
peakfit('GrIr1740CD2longanneal',kalifit.p1,kalifit.p2,'g',2,740,1)
figure
peakfit('GrIr1740CD2_2504',kalifit.p1,kalifit.p2,'r',2,750,1)

figure
peakfit('050516IrSorenDdose1hour2',kalifit.p1,kalifit.p2,'b',2,500,1)

%% Integrate graph. Return area


end
function  TPDplot(filename,kal1,kal2) 
a = filename;
b = '.csv';
TPDdata = importdata([a b],',',36);

Timems = TPDdata.data(:,1)*10^-3;
D2counts = TPDdata.data(:,2);
temp = TPDdata.data(:,3)*kal1+kal2;



plot(Timems,D2counts,'b')
yyaxis left
ylabel('Counts')
hold on
yyaxis right
ylabel('Temp [K]')
xlabel('Time [s]')
plot(Timems,temp,'g','LineWidth',2)
set(gca,'FontSize',14,'FontWeight','bold')

legend('Counts','Temp [K]','Location','NorthWest')


print(filename,'-dpng','-r200')
end
function TPDplottemp(filename,kal1,kal2,color,cut)
a = filename;
b = '.csv';
TPDdata = importdata([a b],',',36);

Times = TPDdata.data(:,1)*10^-3;
D2counts = TPDdata.data(:,2);
temp = TPDdata.data(:,3)*kal1+kal2;

cutoffmax = find(temp == max(temp));
cutoffmin = find(temp == min(temp));

if length(cutoffmin)>1
cutoffmin = cutoffmin(length(cutoffmin));
end

cuttime = Times(find(Times >= cut));

background = mean(D2counts(find(Times == cuttime(1)):cutoffmin));
relevantdata = TPDdata.data(cutoffmin:cutoffmax,2);

D2counts = relevantdata'-background*ones(1,length(relevantdata));
temp = TPDdata.data(cutoffmin:cutoffmax,3);
    
plot(temp,D2counts,color,'LineWidth',1)
ylabel('Counts')
xlabel('Temp [K]')
set(gca,'FontSize',14,'FontWeight','bold')

print([filename 'temp'],'-dpng','-r200')
end
function peakfit(filename,kal1,kal2,color,cut,peakguess,peaknumber)
a = filename;
b = '.csv';
TPDdata = importdata([a b],',',36);

Times = TPDdata.data(:,1)*10^-3;
D2counts = TPDdata.data(:,2);
temp = TPDdata.data(:,3)*kal1+kal2;


cutoffmax = find(temp == max(temp));
cutoffmin = find(temp == min(temp));

if length(cutoffmin)>1
cutoffmin = cutoffmin(length(cutoffmin));
end

cuttime = Times(find(Times >= cut));

background = mean(D2counts(find(Times == cuttime(1)):cutoffmin));
relevantdata = TPDdata.data(cutoffmin:cutoffmax,2);

searchcenter = find(abs(temp-peakguess) == min(abs(temp-peakguess)));

if length(D2counts) < searchcenter+235
Peak = temp(find( D2counts == max(D2counts(searchcenter-235:length(D2counts)))));
elseif searchcenter-235 < 0 
Peak = temp(find( D2counts == max(D2counts(1:searchcenter+235))));
else
Peak = temp(find( D2counts == max(D2counts(searchcenter-235:searchcenter+235))));
end

fitdata = find(abs(temp-Peak(1)) == min(abs(temp-Peak(1))));

fittemp = temp(fitdata(1)-235:fitdata(1)+235);
fitD2counts = D2counts(fitdata(1)-235:fitdata(1)+235);

plot(fittemp,fitD2counts,'xk')
hold on
gaussfit = fit(fittemp,fitD2counts,'gauss3');
h = plot(gaussfit,'r');
set(h, 'LineWidth',2)

xdata = linspace(Peak(1)-40,Peak(1)+40,1000);
[peakvalue peaktemp] = findpeaks(gaussfit(xdata),'NPeaks',peaknumber,'SortStr','descend');
peaktemperature = xdata(peaktemp)

plot(peaktemperature(peaknumber),peakvalue(peaknumber),'og','MarkerSize',15,'LineWidth',2)
ylabel('Counts')
xlabel('Temperature [K]')
set(gca,'FontSize',14,'FontWeight','bold')

print([filename 'peak'],'-dpng','-r200')
end