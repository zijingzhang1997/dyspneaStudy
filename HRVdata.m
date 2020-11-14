 dataPath=['C:\Sleep test\dyspnea\data\10_6\'];
fileName1_1='Case1Routine1';toff1_1=[10:175]';
[hrvFeature1_1,pk1_1] = HRV_ecg(dataPath,fileName1_1,toff1_1);
fileName1_5a='Case1Routine5a';toff1_5a=[10:175]';
[hrvFeature1_5a,pk1_5a] = HRV_ecg(dataPath,fileName1_5a,toff1_5a);

fileName2_1='Case2Routine1';toff2_1=[10:175]';
[hrvFeature2_1,pk2_1] = HRV_ecg(dataPath,fileName2_1,toff2_1);
fileName2_5a='Case2Routine5a';toff2_5a=[10:175]';
[hrvFeature2_5a,pk2_5a] = HRV_ecg(dataPath,fileName2_5a,toff2_5a);

fileName3_1='Case3Routine1';toff3_1=[10:175]';
[hrvFeature3_1,pk3_1] = HRV_ecg(dataPath,fileName3_1,toff3_1);
fileName3_5a='Case3Routine5a';toff3_5a=[10:175]';
[hrvFeature3_5a,pk3_5a] = HRV_ecg(dataPath,fileName3_5a,toff3_5a);


fileName1_5b='Case1Routine5b';toff1_5b=[10:175]';
[hrvFeature1_5b,pk1_5b] = HRV_ecg(dataPath,fileName1_5b,toff1_5b);
fileName2_5b='Case1Routine5b';toff2_5b=[10:175]';
[hrvFeature2_5b,pk2_5b] = HRV_ecg(dataPath,fileName2_5b,toff2_5b);
fileName3_5b='Case1Routine5b';toff3_5b=[10:175]';
[hrvFeature3_5b,pk3_5b] = HRV_ecg(dataPath,fileName3_5b,toff3_5b);
fileName2_3='Case2Routine3';toff2_3=[10:175]';
[hrvFeature2_3,pk2_3] = HRV_ecg(dataPath,fileName2_3,toff2_3);
fileName3_3='Case3Routine3';toff3_3=[10:175]';
[hrvFeature3_3,pk3_3] = HRV_ecg(dataPath,fileName3_3,toff3_3);


function [hrvFeature,pk] = HRV_ecg(dataPath,fileName,toff)
%dataPath=['C:\Sleep test\dyspnea\data\10_6\'];

fs=10e3;
fsDS=500;
%toff=[10:175]';
%toff=[10:280]';

fsbio=1e3;
fsDSbio=500;
filePathName = [dataPath,fileName,'.tdms'];
filePathName_m = [dataPath,fileName,'.mat'];
if ~exist(filePathName_m,'file')
   convertTDMS(true,filePathName);
end
load(filePathName_m);

PPG=ConvertedData.Data.MeasuredData(22).Data;  %19 biopac ch1
ECG=ConvertedData.Data.MeasuredData(19).Data;

ds_ppg=resample(PPG,fsDSbio,fsbio);
ds_ecg=resample(ECG,fsDSbio,fsbio);
t = ((0:(length(ds_ecg)-1))/fsDSbio)';
ds_ppgOff=ds_ppg((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);
ds_ecgOff=ds_ecg((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);


tOff=((0:(length(ds_ecgOff)-1))/fsDSbio)';



opts5.filtType = 'LpHp';
opts5.f3db = 4; opts5.fpLP = 20; opts5.fstLP = 25;
ds_ecgfilt = filterLpHp(ds_ecgOff,fsDS,opts5); 
%ds_ecgfilt=ds_ecgOff;
ds_ppgfilt=ds_ppgOff;


[ds_ppgfilt_norm,PS] = mapminmax(ds_ppgfilt');
ds_ppgfilt_norm=ds_ppgfilt_norm';
[ds_ecgfilt_norm,PS] = mapminmax(ds_ecgfilt');
ds_ecgfilt_norm=ds_ecgfilt_norm';

sz=10;
figure()
subplot(2,1,1);
plot(tOff,ds_ppgfilt,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 20])
titleName=[fileName,'  PPG']
title(titleName,'FontSize',sz)


subplot(2,1,2);
plot(tOff,ds_ecgfilt,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 20])
title('ECG','FontSize',sz)


set(gcf,'Position',[100,200,500,350]);
figName = [dataPath,'fig/',fileName,'Ecg_ppg'];
%  saveas(gca,figName,'tif')
% saveas(gca,figName,'fig')




%% HR by ECG
opts4.tWinBR = 15; % Window on which br is estimated
opts4.tWin = 4; % Window for peak detection moving average
opts4.minInterceptDist = 0.15; 
opts4.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts4.calibMinPkRatio = 0.3;

opts5.filtType = 'LpHp';
opts5.f3db = 4; opts5.fpLP = 20; opts5.fstLP = 25;
ecgFilt = filterLpHp(ds_ecgfilt_norm,fsDS,opts5);
plot(ecgFilt)
opts4.tWinHR = 4; % Same for ECG and NCS
opts4.tWin = 0.5;
opts4.minInterceptDist = 0.05;
opts4.minPkHt = 0.05;
%  [hrBio, pkMaxBio] = ecgHR(ds_ecgfilt_norm,fsDS,opts4);
%  figure()
%  plot(tOff,hrBio);

[hrvFeature,pk] = hrvFeatureEstEcg(ds_ecgfilt_norm,fsDS,opts4);
% savePath='C:\Sleep test\dyspnea\data\10_6\mfile\HRVfeature\';
% save([savePath,fileName,'.mat'],'hrvFeature');
end