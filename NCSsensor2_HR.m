dataPath=['F:\Sleep test\Data\330\'];
fileName='2.0P1';

fs=10e3;
fsDS=500;
toff=[10:230]';
filePathName = [dataPath,fileName,'.tdms'];
%convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);
amp1=ConvertedData.Data.MeasuredData(3).Data;
ph1=ConvertedData.Data.MeasuredData(5).Data;
amp2=ConvertedData.Data.MeasuredData(4).Data;
ph2=ConvertedData.Data.MeasuredData(6).Data;


ampds1=resample(amp1,fsDS,fs);
phds1=resample(ph1,fsDS,fs);
t = ((0:(length(ampds1)-1))/fsDS)';
ampDsOff1=ampds1((toff(1)*fsDS):toff(size(toff))*fsDS);
phDsOff1=phds1((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOff1)-1))/fsDS)';
opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.7; opts1.fpLP = 1.2; opts1.fstLP = 1.5;
ampfilt1 = filterLpHp(ampDsOff1,fsDS,opts1); % th amp
phfilt1 = filterLpHp(phDsOff1,fsDS,opts1); 

ampds2=resample(amp2,fsDS,fs);
phds2=resample(ph2,fsDS,fs);
t = ((0:(length(ampds2)-1))/fsDS)';
ampDsOff2=ampds2((toff(1)*fsDS):toff(size(toff))*fsDS);
phDsOff2=phds2((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOff2)-1))/fsDS)';

opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.7; opts1.fpLP = 1.2; opts1.fstLP = 1.5;
ampfilt2 = filterLpHp(ampDsOff2,fsDS,opts1); % th amp
phfilt2 = filterLpHp(phDsOff2,fsDS,opts1); 

yampDsOff1=fft(ampfilt1);
f = (0:length(yampDsOff1)-1)*fsDS/length(yampDsOff1)';
yam1=abs(yampDsOff1);
f=f(1,1:1000);
yam1=yam1(1:1000,1)';
yam1max=max(yam1);
yam1min=min(yam1);
yamnorm1=zeros(1,1000);
yamnorm1=(yam1-yam1min)/(yam1max-yam1min);
figure
plot(f,yamnorm1,'LineWidth',0.5,'color','blue');
xlabel('frequency/Hz')

title('amp th ')

yampDsOff2=fft(ampfilt2);
f = (0:length(yampDsOff2)-1)*fsDS/length(yampDsOff2)';
yam2=abs(yampDsOff2);
f=f(1,1:1000);
yam2=yam2(1:1000,1)';
yam2max=max(yam2);
yam2min=min(yam2);
yamnorm2=zeros(1,1000);
yamnorm2=(yam2-yam2min)/(yam2max-yam2min);
figure
plot(f,yamnorm2,'LineWidth',0.5,'color','blue');
xlabel('frequency/Hz')

title('amp ab ')


figure

subplot(2,1,1);
plot(tOff,ampfilt1,'LineWidth',0.5,'color','red');
xlabel('t/s')

 title('Amp TH filtered HR')

subplot(2,1,2);
plot(tOff,ampfilt2,'LineWidth',0.5,'color','green');
xlabel('t/s')
title('Amp AB filtered HR')