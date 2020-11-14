dataPath=['F:\Sleep test\Data\301positionbias_th\'];
fileName='0.900absideL';
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
opts1.f3db = 0.05; opts1.fpLP = 1; opts1.fstLP = 1.5;
ampfilt = filterLpHp(ampDsOff,fsDS,opts1); 
phfilt = filterLpHp(phDsOff,fsDS,opts1); 





% 
% figure
% 
% subplot(2,1,1);
% plot(tOff,ampfilt,'LineWidth',0.5,'color','red');
% xlabel('t/s')
%  title('amp 1.8GHz supine ')
% 
% subplot(2,1,2);
% plot(tOff,phfilt,'LineWidth',0.5,'color','green');
% xlabel('t/s')
% title('ph 1.8GHz supine ')


 opts3.tWinBR = 16.5; % Window on which br is estimated
opts3.tWin = 4; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.3;

[br,RespPk] = brEst(ampfilt,fsDS,opts3);
pkMax = RespPk.idx(RespPk.ind == 1);
pkMin = RespPk.idx(RespPk.ind == 0);
pkMax = pkMax(RespPk.idxValidPk);
pkMin = pkMin(RespPk.idxValidPk);
tpk=tOff(pkMax);
pkValue=ampfilt(pkMax)-ampfilt(pkMin);
couple=pkValue./ampDsOff(pkMax);
couple=100*couple;

figure()
plot(tpk,couple)
ylabel('coupling strength(%)')

% ampfilt1max=max(ampfilt1);
% ampfilt1min=min(ampfilt1);
% ampfilt1norm=zeros(1,length(ampfilt1));
% ampfilt1norm=(ampfilt1-ampfilt1min)/(ampfilt1max-ampfilt1min);
% phfilt1max=max(phfilt1);
% phfilt1min=min(phfilt1);
% phfilt1norm=zeros(1,length(phfilt1));
% phfilt1norm=(phfilt1-phfilt1min)/(phfilt1max-phfilt1min);
% ampfilt2max=max(ampfilt2);
% ampfilt2min=min(ampfilt2);
% ampfilt2norm=zeros(1,length(ampfilt2));
% ampfilt2norm=(ampfilt2-ampfilt2min)/(ampfilt2max-ampfilt2min);
% phfilt2max=max(phfilt2);
% phfilt2min=min(phfilt2);
% phfilt2norm=zeros(1,length(phfilt2));
% phfilt2norm=(phfilt2-phfilt2min)/(phfilt2max-phfilt2min);
% 
% figure()
% plot(tOff,ampfilt1norm,'LineWidth',0.5,'color','red');
% 
% hold on
% plot(tOff,-ampfilt2norm,'LineWidth',0.5,'color','green');
% xlabel('t/s')
% ylabel('normalized Amplitude')
% legend('amp TH','amp AB')


