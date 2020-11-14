dataPath=['F:\Sleep test\Data\226positionbias_ab\'];
fileName1='0.900';
fileName2='0.904';
fileName3='0.908';
fs=25e3;
fsDS=500;
toff=[10:170]';
%filePathName = [dataPath,fileName1,'.tdms'];
%convertTDMS(true,filePathName);

load([dataPath,fileName1,'.mat']);
amp1=ConvertedData.Data.MeasuredData(3).Data;
load([dataPath,fileName2,'.mat']);
amp2=ConvertedData.Data.MeasuredData(3).Data;
load([dataPath,fileName3,'.mat']);
amp3=ConvertedData.Data.MeasuredData(3).Data;

ampds1=resample(amp1,fsDS,fs);  
t = ((0:(length(ampds1)-1))/fsDS)';
ampDsOff1=ampds1((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOff1)-1))/fsDS)';
ampds2=resample(amp2,fsDS,fs);  
ampDsOff2=ampds2((toff(1)*fsDS):toff(size(toff))*fsDS);
ampds3=resample(amp3,fsDS,fs);  
ampDsOff3=ampds3((toff(1)*fsDS):toff(size(toff))*fsDS);

opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 1; opts1.fstLP = 1.5;
ampfilt1 = filterLpHp(ampDsOff1,fsDS,opts1); 
ampfilt2 = filterLpHp(ampDsOff2,fsDS,opts1); 
ampfilt3 = filterLpHp(ampDsOff3,fsDS,opts1); 
sz=10;
[ampfilt1_norm,PS] = mapminmax(ampfilt1');
ampfilt1_norm=ampfilt1_norm';
[ampfilt2_norm,PS] = mapminmax(ampfilt2');
ampfilt2_norm=ampfilt2_norm';
[ampfilt3_norm,PS] = mapminmax(ampfilt3');
ampfilt3_norm=ampfilt3_norm';

 opts3.tWinBR = 16.5; % Window on which br is estimated
opts3.tWin = 4; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.4;
[br1,RespPk1] = brEst(ampfilt1,fsDS,opts3);
[br2,RespPk1] = brEst(ampfilt2,fsDS,opts3);
[br3,RespPk1] = brEst(ampfilt3,fsDS,opts3);




figure()
subplot(3,1,1)
plot(tOff,ampfilt1_norm,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 160])
title('d=0 (inch) 2.4GHz','FontSize',sz)

subplot(3,1,2)
plot(tOff,ampfilt2_norm,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 160])
title('d=4 (inch) 2.4GHz','FontSize',sz)

subplot(3,1,3)
plot(tOff,ampfilt3_norm,'LineWidth',0.5,'color','blue');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 160])
title('d=8 (inch) 2.4GHz','FontSize',sz)




figure()
plot(tOff,br1,'LineWidth',1,'color','gree');
xlabel('time(s)','FontSize',sz)
ylabel('BR(BPM)','FontSize',sz)
hold on 
plot(tOff,br2,'LineWidth',1,'color','red');
hold on 
plot(tOff,br3,'LineWidth',1,'color','blue');

xlabel('time (s)','FontSize',sz)
ylabel('BR (BPM)','FontSize',sz)
legend('d=0','d=4','d=8','FontSize',sz)
% pkMax = RespPk.idx(RespPk.ind == 1);
% pkMin = RespPk.idx(RespPk.ind == 0);
% pkMax = pkMax(RespPk.idxValidPk);
% pkMin = pkMin(RespPk.idxValidPk);
% tpk=tOff(pkMax);
% pkValue=ampfilt(pkMax)-ampfilt(pkMin);
% couple=pkValue./ampDsOff(pkMax);
% couple=100*couple;
% f3=figure()
% plot(tpk,couple)
% xlabel('time(s)','FontSize',sz)
% ylabel('coupling strength(%)','FontSize',sz)




