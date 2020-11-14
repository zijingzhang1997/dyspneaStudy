dataPath=['F:\Sleep test\Data\317NCS2bio2\'];
fileName='1.8';
% fs=25e3;
fs=10e3;
fsDS=500;
toff=[10:230]';
sz=10;
filePathName = [dataPath,fileName,'.tdms'];
 convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);
amp=ConvertedData.Data.MeasuredData(3).Data;

ph=ConvertedData.Data.MeasuredData(3).Data;


ampds=resample(amp,fsDS,fs);
phds=resample(ph,fsDS,fs);
t = ((0:(length(ampds)-1))/fsDS)';
ampDsOff=ampds((toff(1)*fsDS):toff(size(toff))*fsDS);
phDsOff=phds((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOff)-1))/fsDS)';

opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.8; opts1.fpLP = 2; opts1.fstLP = 3;
ampfilt = filterLpHp(ampDsOff,fsDS,opts1); % th amp
phfilt = filterLpHp(phDsOff,fsDS,opts1); 

yampDsOff=fft(ampfilt);
f = (0:length(yampDsOff)-1)*fsDS/length(yampDsOff)';
yam=abs(yampDsOff);
f=f(1,1:1000);
yam=yam(1:1000,1)';
yammax=max(yam);
yammin=min(yam);
yamnorm=zeros(1,1000);
yamnorm=(yam-yammin)/(yammax-yammin);


% yphDsOff=fft(phfilt);
% f = (0:length(yphDsOff)-1)*fsDS/length(yphDsOff)';
% yph=abs(yphDsOff);
% f=f(1,1:100);
% yph=yph(1:100,1)';
% yphmax=max(yph);
% yphmin=min(yph);
% yphnorm=zeros(1,100);
% yphnorm=(yph-yphmin)/(yphmax-yphmin);
% figure
% plot(f,yphnorm,'LineWidth',0.5,'color','blue');
% xlabel('frequency/Hz')
% title('ph ');

[ampfilt,PS] = mapminmax(ampfilt');
figure
subplot(1,2,2)
plot(f,yamnorm,'LineWidth',0.8,'color','red');
xlabel('Frequency (Hz)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 4]);

subplot(1,2,1)
plot(tOff,ampfilt','LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
title('1.8GHz','FontSize',sz)
xlim([0 220]);

% 
% subplot(2,2,3);
% plot(f,y,'LineWidth',0.5,'color','blue');
% title('amp ');
% 
% subplot(2,2,4);
% plot(tOff,phDsOff,'LineWidth',0.5,'color','blue');
% title('ph ');


% opts3.tWinBR = 15; % Window on which br is estimated
% opts3.tWin = 8; % Window for peak detection moving average
% opts3.minInterceptDist = 0.15; 
% opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
% opts3.calibMinPkRatio = 0.6;
% 
% [brNcs,ncsRespPk] = brEst(ampfilt,fsDS,opts3);
% pkMaxNcs = ncsRespPk.idx(ncsRespPk.ind == 1);
% pkMinNcs = ncsRespPk.idx(ncsRespPk.ind == 0);
% pkMaxNcs = pkMaxNcs(ncsRespPk.idxValidPk);
% pkMinNcs = pkMinNcs(ncsRespPk.idxValidPk);
% 
% 
% [brNcs,ncsRespPk] = brEst(phfilt,fsDS,opts3);
% pkMaxNcs = ncsRespPk.idx(ncsRespPk.ind == 1);
% pkMinNcs = ncsRespPk.idx(ncsRespPk.ind == 0);
% pkMaxNcs = pkMaxNcs(ncsRespPk.idxValidPk);
% pkMinNcs = pkMinNcs(ncsRespPk.idxValidPk);

