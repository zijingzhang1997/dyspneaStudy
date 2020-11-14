dataPath=['C:\Sleep test\Data\603\'];
fileName='250';
fs=10e3;
fsDS=500;
toff=[30:120]';
filePathName = [dataPath,fileName,'.tdms'];
convertTDMS(true,filePathName);


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
opts1.f3db = 0.05; opts1.fpLP = 0.8; opts1.fstLP = 1;
ampfilt1 = filterLpHp(ampDsOff1,fsDS,opts1); % th amp
phfilt1 = filterLpHp(phDsOff1,fsDS,opts1); 

ampds2=resample(amp2,fsDS,fs);
phds2=resample(ph2,fsDS,fs);
t = ((0:(length(ampds2)-1))/fsDS)';
ampDsOff2=ampds2((toff(1)*fsDS):toff(size(toff))*fsDS);
phDsOff2=phds2((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOff2)-1))/fsDS)';

opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 0.8 ;opts1.fstLP = 1;
ampfilt2 = filterLpHp(ampDsOff2,fsDS,opts1); % th amp
phfilt2 = filterLpHp(phDsOff2,fsDS,opts1); 


figure

subplot(4,1,1);
plot(tOff,ampfilt1,'LineWidth',0.5,'color','red');
xlabel('t/s')

 title('Amp TH')

subplot(4,1,2);
plot(tOff,phfilt1,'LineWidth',0.5,'color','green');
xlabel('t/s')
title('Ph TH')

subplot(4,1,3);
plot(tOff,ampfilt2,'LineWidth',0.5,'color','red');
xlabel('t/s')
 title('Amp AB')

subplot(4,1,4);
plot(tOff,phfilt2,'LineWidth',0.5,'color','green');
xlabel('t/s')
title('Ph AB')

ampfilt1max=max(ampfilt1);
ampfilt1min=min(ampfilt1);
ampfilt1norm=zeros(1,length(ampfilt1));
ampfilt1norm=(ampfilt1-ampfilt1min)/(ampfilt1max-ampfilt1min);
phfilt1max=max(phfilt1);
phfilt1min=min(phfilt1);
phfilt1norm=zeros(1,length(phfilt1));
phfilt1norm=(phfilt1-phfilt1min)/(phfilt1max-phfilt1min);
ampfilt2max=max(ampfilt2);
ampfilt2min=min(ampfilt2);
ampfilt2norm=zeros(1,length(ampfilt2));
ampfilt2norm=(ampfilt2-ampfilt2min)/(ampfilt2max-ampfilt2min);
phfilt2max=max(phfilt2);
phfilt2min=min(phfilt2);
phfilt2norm=zeros(1,length(phfilt2));
phfilt2norm=(phfilt2-phfilt2min)/(phfilt2max-phfilt2min);

figure()
plot(tOff,ampfilt1norm,'LineWidth',0.5,'color','red');

hold on
plot(tOff,ampfilt2norm,'LineWidth',0.5,'color','green');
xlabel('t/s')
ylabel('normalized Amplitude')
legend('amp TH','amp AB')

