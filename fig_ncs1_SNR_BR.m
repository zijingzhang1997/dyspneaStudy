dataPath=['F:\Sleep test\Data\226positionbias_ab\'];
fileName='1.800';
fs=25e3;
fsDS=500;
toff=[10:55]';
filePathName = [dataPath,fileName,'.tdms'];
%convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);
amp=ConvertedData.Data.MeasuredData(3).Data;
ph=ConvertedData.Data.MeasuredData(4).Data;


ampds=resample(amp,fsDS,fs);  
phds=resample(ph,fsDS,fs);
t = ((0:(length(ampds)-1))/fsDS)';
ampDsOff=ampds((toff(1)*fsDS):toff(size(toff))*fsDS);
phDsOff=phds((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOff)-1))/fsDS)';

opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 1; opts1.fstLP = 1.5;
ampfilt = filterLpHp(ampDsOff,fsDS,opts1); 
phfilt = filterLpHp(phDsOff,fsDS,opts1); 
sz=11;
% [ampfilt_norm,PS] = mapminmax(ampfilt');
% ampfilt_norm=ampfilt_norm';

f1=figure()
plot(tOff,ampfilt,'LineWidth',0.5,'color','red');
xlabel('time(s)','FontSize',sz)
ylabel('Amp(a.u.)','FontSize',sz)
xlim([0 160])
% ylim([-3 3])
title('x=8(inch)','FontSize',sz)

yampDsOff=fft(ampfilt);
f = (0:length(yampDsOff)-1)*fsDS/length(yampDsOff)';
yam=abs(yampDsOff);
f=f(1,1:60);
yam=yam(1:60,1)';
yammax=max(yam);
yammin=min(yam);
yamnorm=zeros(1,60);
yamnorm=(yam-yammin)/(yammax-yammin);
figure
plot(f,yamnorm,'LineWidth',0.8,'color','blue');
xlabel('Frequency (Hz)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 3]);


snr(ampfilt,fsDS);
yampDsOff=fft(ampfilt);
f = (0:length(yampDsOff)-1)*fsDS/length(yampDsOff)';
yam=abs(yampDsOff);
f=f(1,1:60);
yam=yam(1:60,1)';
yammax=max(yam);
yammin=min(yam);
yamnorm=zeros(1,60);
yamnorm=(yam-yammin)/(yammax-yammin);
yamnorm=yam;
figure
plot(f,yamnorm,'LineWidth',0.8,'color','blue');
xlabel('Frequency (Hz)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 3]);
yamnormdb=mag2db(yamnorm);
figure
plot(f,yamnormdb,'LineWidth',0.8,'color','blue');
xlabel('Frequency (Hz)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 3]);


