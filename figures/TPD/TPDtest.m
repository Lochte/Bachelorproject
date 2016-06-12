function TPDtest()
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


%% 
figure
TPDplotarhennius('GrIr1740CD2_2504',kalifit.p1,kalifit.p2,'rx',2,753.4)

TPDplotarhennius('GrIr1740CD2longanneal',kalifit.p1,kalifit.p2,'gx',2,739.1)

TPDplotarhennius('050516IrSorenD2dose1hour2',kalifit.p1,kalifit.p2,'bx',2,724.8)

TPDplotarhennius('050516IrSorenD2dose1hour',kalifit.p1,kalifit.p2,'gx',2,698.4)

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
function TPDplotarhennius(filename,kal1,kal2,color,cut,peaktemp)
a = filename;
b = '.csv';
TPDdata = importdata([a b],',',36);

Times = TPDdata.data(:,1)*10^-3;
D2counts = TPDdata.data(:,2);
temp = TPDdata.data(:,3)*kal1+kal2;

Kb = 8.6173324*10^-5; %eV K^-1
Kb2 = 1.38064852*10^-23;

cutoffmax = find(temp >= peaktemp-10);
cutoffmin = find(temp >= min(temp)+250);

if length(cutoffmax)>1
cutoffmax = cutoffmax(1);
end

if length(cutoffmin)>1
cutoffmin = cutoffmin(1);
end

cuttime = Times(find(Times >= cut));

background = mean(D2counts(find(Times == cuttime(1)):cutoffmin));
relevantdata = TPDdata.data(cutoffmin:cutoffmax,2);

D2counts = relevantdata'-background*ones(1,length(relevantdata));
temp = TPDdata.data(cutoffmin:cutoffmax,3);

totalH = sum(D2counts);

normD2counts = cumsum(D2counts/totalH);

theta = 1-normD2counts;

figure
plot(temp,theta)
hold on
plot(temp(1:length(-diff(theta))),-diff(theta))




figure
plot(power(temp*Kb,-1),log(D2counts/trapz(temp,D2counts)),color,'LineWidth',1)
arrhfit = fit(power(temp*Kb,-1),log(D2counts/trapz(temp,D2counts))','poly1');
hold on
plot(arrhfit)

figure 
plot(temp,D2counts)
Edesorp = -arrhfit.p1
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