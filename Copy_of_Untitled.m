dataPath=['F:\Sleep test\Data\301\'];
fileName='2.404th';
fs=25e3;
fsDS=500;
toff=[10:170]';
filePathName = [dataPath,fileName,'.tdms'];
 convertTDMS(true,filePathName);


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
opts1.f3db = 0.2; opts1.fpLP = 5; opts1.fstLP = 7;
ampfilt = filterLpHp(ampDsOff,fsDS,opts1); % th amp
phfilt = filterLpHp(phDsOff,fsDS,opts1); 


yampDsOff=fft(ampfilt);
f = (0:length(yampDsOff)-1)*fsDS/length(yampDsOff)';
yam=abs(yampDsOff);
f=f(1,1:100);
yam=yam(1:100,1)';
yammax=max(yam);
yammin=min(yam);
yamnorm=zeros(1,100);
yamnorm=(yam-yammin)/(yammax-yammin);
% figure
% plot(f,yamnorm,'LineWidth',0.5,'color','blue');
% xlabel('frequency/Hz')

% title('amp ')

yphDsOff=fft(phfilt);
f = (0:length(yphDsOff)-1)*fsDS/length(yphDsOff)';
yph=abs(yphDsOff);
f=f(1,1:100);
yph=yph(1:100,1)';
yphmax=max(yph);
yphmin=min(yph);
yphnorm=zeros(1,100);
yphnorm=(yph-yphmin)/(yphmax-yphmin);
% figure
% plot(f,yphnorm,'LineWidth',0.5,'color','blue');
% xlabel('frequency/Hz')
% title('ph ');

figure

subplot(2,1,1);
plot(tOff,ampfilt,'LineWidth',0.5,'color','red');
xlabel('t/s')
 title('amp  2.4GHz chest x=8inch h=7inch')

subplot(2,1,2);
plot(tOff,phfilt,'LineWidth',0.5,'color','green');
xlabel('t/s')
title('ph  2.4GHz chest x=8inch h=7inch')
% 
% subplot(2,2,3);
% plot(f,y,'LineWidth',0.5,'color','blue');
% title('amp ');
% 
% subplot(2,2,4);
% plot(tOff,phDsOff,'LineWidth',0.5,'color','blue');
% title('ph ');


opts3.tWinBR = 15; % Window on which br is estimated
opts3.tWin = 8; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.6;

[brNcs,ncsRespPk] = brEst(ampfilt,fsDS,opts3);
pkMaxNcs = ncsRespPk.idx(ncsRespPk.ind == 1);
pkMinNcs = ncsRespPk.idx(ncsRespPk.ind == 0);
pkMaxNcs = pkMaxNcs(ncsRespPk.idxValidPk);
pkMinNcs = pkMinNcs(ncsRespPk.idxValidPk);


[brNcs,ncsRespPk] = brEst(phfilt,fsDS,opts3);
pkMaxNcs = ncsRespPk.idx(ncsRespPk.ind == 1);
pkMinNcs = ncsRespPk.idx(ncsRespPk.ind == 0);
pkMaxNcs = pkMaxNcs(ncsRespPk.idxValidPk);
pkMinNcs = pkMinNcs(ncsRespPk.idxValidPk);

