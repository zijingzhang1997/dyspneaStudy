dataPath=['F:\Sleep test\Data\317\'];
fileName='2.4SIDER_0317_134021Routine1'
fs=10e3;
fsDS=500;
toff=[20:230]';
fsbio=1e3;
fsDSbio=500;
filePathName = [dataPath,fileName,'.tdms'];
convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);
ampTh=ConvertedData.Data.MeasuredData(8).Data;
bioTh=ConvertedData.Data.MeasuredData(20).Data;
ampAb=ConvertedData.Data.MeasuredData(3).Data;
bioAb=ConvertedData.Data.MeasuredData(21).Data;


ampdsTh=resample(ampTh,fsDS,fs);
ampdsAb=resample(ampAb,fsDS,fs);
t = ((0:(length(ampdsTh)-1))/fsDS)';
ampDsOffTh=ampdsTh((toff(1)*fsDS):toff(size(toff))*fsDS);
ampDsOffAb=ampdsAb((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOffTh)-1))/fsDS)';


dsbioTh=resample(bioTh,fsDSbio,fsbio);
dsbioAb=resample(bioAb,fsDSbio,fsbio);
t = ((0:(length(dsbioTh)-1))/fsDSbio)';
dsbioThOff=dsbioTh((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);
dsbioAbOff=dsbioAb((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);
tOff=((0:(length(dsbioThOff)-1))/fsDSbio)';

% opts1.filtType = 'LpHp'; opts1.orderHP = 5;
% opts1.f3db = 0.05; opts1.fpLP = 5; opts1.fstLP = 7;



opts2.filtType = 'LpHp'; opts2.orderHP = 5;
opts2.f3db = 0.7; opts2.fpLP = 2; opts2.fstLP = 2.5;
ampfiltHRTh = filterLpHp(ampDsOffTh,fsDS,opts2); 
ampfiltHRAb = filterLpHp(ampDsOffAb,fsDS,opts2); 

figure

subplot(2,1,1);
plot(tOff,ampfiltHRTh,'LineWidth',0.5,'color','red');
xlabel('t/s')

 title('Amp NCS TH filtered heartbeat')

subplot(2,1,2);
plot(tOff,ampfiltHRAb,'LineWidth',0.5,'color','green');
xlabel('t/s')
title('Amp NCS Ab filtered heartbeat')



%heart beat rate filter
yampDsOff=fft(ampfiltHRTh);
f = (0:length(yampDsOff)-1)*fsDS/length(yampDsOff)';
yam=abs(yampDsOff);
f=f(1,1:1000);
yam=yam(1:1000,1)';
yammax=max(yam);
yammin=min(yam);
yamnorm=zeros(1,1000);
yamnorm=(yam-yammin)/(yammax-yammin);
figure
plot(f,yamnorm,'LineWidth',0.5,'color','blue');
xlabel('frequency/Hz')

title('amp Th  frequency spectrum')

yampAbDsOff=fft(ampfiltHRAb);
f = (0:length(yampAbDsOff)-1)*fsDS/length(yampAbDsOff)';
yamAb=abs(yampAbDsOff);
f=f(1,1:1000);
yamAb=yamAb(1:1000,1)';
yamAbmax=max(yamAb);
yamAbmin=min(yamAb);
yamAbnorm=zeros(1,1000);
yamAbnorm=(yamAb-yamAbmin)/(yamAbmax-yamAbmin);
figure
plot(f,yamAbnorm,'LineWidth',0.5,'color','blue');
xlabel('frequency/Hz')

title('amp Ab  frequency spectrum')