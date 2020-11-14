dataPath=['F:\Sleep test\Data\326\'];
fileName='0.9P2';
fs=10e3;
fsDS=500;
toff=[20:240]';
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


slopeBNTA1 = diff(ampfilt1norm).*fsDS; 
slopeBNTA1 = [slopeBNTA1(1,:); slopeBNTA1]; % Same row size as t

% Rescaling slope, for cases when it shoots up a lot
slopeBNTA1 = tanh(slopeBNTA1);
slopeBNTA2 = diff(ampfilt2norm).*fsDS; 
slopeBNTA2 = [slopeBNTA2(1,:); slopeBNTA2]; % Same row size as t

% Rescaling slope, for cases when it shoots up a lot
slopeBNTA2 = tanh(slopeBNTA2);

numFeatIV = 3; 
featDescriptIV = {'Avg SP','StdDev SP','Correlation'};
% Features: [slope product, Avg slope product, Std Dev slope product, correlation] 
% Moving average and standard deviation are taken on a window of length
% tWinFeat, and updated only after tWinSlideFeat
tWinFeat = 3; % Estimate (avg, stddev) features over tWinFeat seconds
tWinCorrFeat = 3; % (s) Window for correlation feat, centered at current index.
tWinSlideFeat = 3; % This is the duration between each window in seconds
numWinFeat = floor((tOff(end)-tWinFeat)/tWinSlideFeat)+1;


slopeProdNcs = (slopeBNTA1.*slopeBNTA2);

movAvgNcsSP = movavg(slopeProdNcs,'Linear',tWinFeat*fsDS); % Moving Avg over past tWinFeat s
movStdNcsSP = movstd(slopeProdNcs,[tWinFeat*fsDS-1,0]); % Moving StdDev over past tWinFeat s





featNcsIV = zeros(numWinFeat,numFeatIV);
tIdxFeatIV = zeros(numWinFeat,1);
for iter = 1:numWinFeat
    % Just keep every ith index of the estimated features, starting from
    % the end of first window index (as moving average and moving standard
    % deviation are estimated with current index and past samples as
    % window.
    idx = (iter)*tWinSlideFeat*fsDS;
    if idx > length(tOff)
        idx = length(tOff);
    end
   

    featNcsIV(iter,1) = movAvgNcsSP(idx);
    featNcsIV(iter,2) = movStdNcsSP(idx);
    %featNcsIV(iter,3) = corrBioNcsTA(idx,2);
    tIdxFeatIV(iter) = idx;
end

tIsoFeat = tOff(tIdxFeatIV);
% -------------------------------------------------------------------------
% Plotting features on previous graph


figure()
plot(tIsoFeat,featNcsIV(:,1));
figure()
plot(tIsoFeat,featNcsIV(:,2));

% -------------------------------------------------------------------------
legend(ax(1),{'Thorax','Abdomen','MA','MStdDev','Corr'},'FontSize',fSz,...
    'Location','southwest','Orientation','Horizontal');
idxCalibIV = ((tIsoFeat >= tCalibIsoVol(1))&(tIsoFeat <= tCalibIsoVol(2)));

% Calibration: as mean feature over the normal breathing period
calibIVFeatBio = mean(featBioIV(idxCalibIV,:)); % Mean of each column
calibIVFeatNcs = mean(featNcsIV(idxCalibIV,:)); 

% Classify as Iso-Volumetric or Paradoxical motion of rib-cage and abdomen
% Using Moving Average as feature of choice, this can be further improved.

if calibIVFeatBio(1)<0
    threshIVBio = 1.05.*calibIVFeatBio;
else
    threshIVBio =  0.*calibIVFeatBio;
end
if calibIVFeatNcs(1)<0
    threshIVNcs = 1.05.*calibIVFeatNcs;
else
    threshIVNcs = 0.*calibIVFeatNcs;
end

indIVBio = featBioIV(:,1) < threshIVBio(1); % Assuming MA is positive otherwise, not using calibration
indIVNcs = featNcsIV(:,1) < threshIVNcs(1);

indIVBioFalse = indIVBio; % Made 0 where there is true IV later
indIVNcsFalse = indIVNcs;

bioIVLevelMA = zeros(numWinFeat,1); % Only MA feature: 0 where below threshold, otherwise MA 
ncsIVlevelMA = zeros(numWinFeat,1); 
bioIVLevelMA(indIVBio) = featBioIV(indIVBio,1);
ncsIVlevelMA(indIVNcs) = featNcsIV(indIVNcs,1);