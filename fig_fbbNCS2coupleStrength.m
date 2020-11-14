f=[500,1e3,5e3,10e3,20e3,50e3,100e3];
f=log10(f);

fileName1='500';fileName2='1k';fileName3='5k';fileName4='10K2';fileName5='20k';fileName6='50k (2)';fileName7='100K2 (2)';
[std_cpTh(1),mean_cpTh(1),std_cpAb(1),mean_cpAb(1),snrTh(1),snrAb(1)]=couple_ratio(fileName1);
[std_cpTh(2),mean_cpTh(2),std_cpAb(2),mean_cpAb(2),snrTh(2),snrAb(2)]=couple_ratio(fileName2);
[std_cpTh(3),mean_cpTh(3),std_cpAb(3),mean_cpAb(3),snrTh(3),snrAb(3)]=couple_ratio(fileName3);
[std_cpTh(4),mean_cpTh(4),std_cpAb(4),mean_cpAb(4),snrTh(4),snrAb(4)]=couple_ratio(fileName4);
[std_cpTh(5),mean_cpTh(5),std_cpAb(5),mean_cpAb(5),snrTh(5),snrAb(5)]=couple_ratio(fileName5);
[std_cpTh(6),mean_cpTh(6),std_cpAb(6),mean_cpAb(6),snrTh(6),snrAb(6)]=couple_ratio(fileName6);
[std_cpTh(7),mean_cpTh(7),std_cpAb(7),mean_cpAb(7),snrTh(7),snrAb(7)]=couple_ratio(fileName7);

figure()
sz=10;

errorbar(f,mean_cpAb,std_cpAb,'LineWidth',1.5);
set(gca,'XScale','log')

xticks([f(1) f(2) f(3) f(4) f(5) f(6) f(7)]);
xticklabels({'5\times10^2','10^3','5\times10^3','10^4','2\times10^4','5\times10^4','10^5'})

xlabel('{\it f_{BB}} (Hz)','FontSize',sz)
ylabel('Signal Quality (bits)','FontSize',sz)

% figure()
% sz=10;
% errorbar(mean_cpTh,std_cpTh,'LineWidth',0.7);
% hold on
% errorbar(mean_cpAb,std_cpAb,'LineWidth',0.7);
% xticks([1 2 3 4 5 6 7]);
% xticklabels({'500','1k','5k','10k','20k','50k','100k'})
% 
% xlabel('Frequency (Hz)','FontSize',sz)
% ylabel('coupling strength(%)','FontSize',sz)
% legend('Th','Ab')
% figure()
% plot(snrTh)
% hold on 
% plot(snrAb)
function [std_cpTh,mean_cpTh,std_cpAb,mean_cpAb,snrTh,snrAb] = couple_ratio(fileName)
dataPath=['C:\Sleep test\Data\519\'];

fs=10e3;
fsDS=500;
toff=[30:120]';
filePathName = [dataPath,fileName,'.tdms'];
%convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);
amp1=ConvertedData.Data.MeasuredData(3).Data;
ph1=ConvertedData.Data.MeasuredData(5).Data;
amp2=ConvertedData.Data.MeasuredData(4).Data;
ph2=ConvertedData.Data.MeasuredData(6).Data;
% amp1=ConvertedData.Data.MeasuredData(8).Data;
% ph1=ConvertedData.Data.MeasuredData(8).Data;
% amp2=ConvertedData.Data.MeasuredData(3).Data;
% ph2=ConvertedData.Data.MeasuredData(3).Data;


ampds1=resample(amp1,fsDS,fs);
phds1=resample(ph1,fsDS,fs);
t = ((0:(length(ampds1)-1))/fsDS)';
ampDsOff1=ampds1((toff(1)*fsDS):toff(size(toff))*fsDS);
phDsOff1=phds1((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOff1)-1))/fsDS)';
opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 1.5; opts1.fstLP = 2;
ampfilt1 = filterLpHp(ampDsOff1,fsDS,opts1); % th amp
phfilt1 = filterLpHp(phDsOff1,fsDS,opts1); 

ampds2=resample(amp2,fsDS,fs);
phds2=resample(ph2,fsDS,fs);
t = ((0:(length(ampds2)-1))/fsDS)';
ampDsOff2=ampds2((toff(1)*fsDS):toff(size(toff))*fsDS);
phDsOff2=phds2((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOff2)-1))/fsDS)';

opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 1.5 ;opts1.fstLP = 2;
ampfilt2 = filterLpHp(ampDsOff2,fsDS,opts1); % th amp
phfilt2 = filterLpHp(phDsOff2,fsDS,opts1); 


% figure()
% 
% subplot(4,1,1);
% plot(tOff,ampfilt1,'LineWidth',0.5,'color','red');
% xlabel('t/s')
% 
%  title('')
%  
% subplot(4,1,2);
% plot(tOff,phfilt1,'LineWidth',0.5,'color','blue');
% xlabel('t/s')
% 
%  title('')
% 
% subplot(4,1,3);
% plot(tOff,ampfilt2,'LineWidth',0.5,'color','red');
% xlabel('t/s')
%  title('')
% 
% subplot(4,1,4);
% plot(tOff,phfilt2,'LineWidth',0.5,'color','blue');
% xlabel('t/s')
%  title('')
 
 %% breath rate and calculate max-min
 
 opts3.tWinBR = 15; % Window on which br is estimated
opts3.tWin = 4; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.3;

[brTh,ThRespPk] = brEst(ampfilt1,fsDS,opts3);
pkMaxTh = ThRespPk.idx(ThRespPk.ind == 1);
pkMinTh = ThRespPk.idx(ThRespPk.ind == 0);
pkMaxTh = pkMaxTh(ThRespPk.idxValidPk);
pkMinTh = pkMinTh(ThRespPk.idxValidPk);
tpkTh=tOff(pkMaxTh);
pkValueTh=ampfilt1(pkMaxTh)-ampfilt1(pkMinTh);
coupleTh=pkValueTh./ampDsOff1(pkMaxTh);
coupleTh=log2(coupleTh)+11;




[brAb,AbRespPk] = brEst(ampfilt2,fsDS,opts3);
pkMaxAb = AbRespPk.idx(AbRespPk.ind == 1);
pkMinAb = AbRespPk.idx(AbRespPk.ind == 0);
pkMaxAb = pkMaxAb(AbRespPk.idxValidPk);
pkMinAb = pkMinAb(AbRespPk.idxValidPk);
tpkAb=tOff(pkMaxAb);
pkValueAb=ampfilt2(pkMaxAb)-ampfilt2(pkMinAb);
coupleAb=pkValueAb./ampDsOff2(pkMaxAb);
coupleAb=log2(coupleAb)+11;

std_cpTh=std(coupleTh);
mean_cpTh=mean(coupleTh);
std_cpAb=std(coupleAb);
mean_cpAb=mean(coupleAb);

snrTh = snr(ampfilt1,fsDS,3);
snrAb = snr(ampfilt2,fsDS,3);
% plot(tpkTh,coupleTh)
% hold on
% plot(tpkAb,coupleAb)
% xlabel('t/s')
% ylabel('coupling strength(%)')
% legend('th','ab')
    
end