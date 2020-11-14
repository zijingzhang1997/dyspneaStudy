
clear all;close all;clc

% This is changed each time, starts from 3, skips 5 (3,8,13...)
fileNum = 108;

% Detection of Isovolumetric breathing:
[caseOrder,routineName,~,~,~] = case7_30Info();
numFiles = length(caseOrder);
numFilePerSubj = 5;
numSubj = numFiles/numFilePerSubj;
dataPath = 'C:\Research\NCS\HumanStudyData\Analysis3_Oct19_final\';

readVar = {'ncs','bio','fs','tNcs','tBio','ncsRespTh','ncsRespAbd'};

caseNum = caseOrder(fileNum);
% fileIdx: 1=Supine, 2=LR, 3=Sitting
fileIdx = mod(fileNum,5);
if fileIdx == 3 
    fprintf('Continue with this file.\n'); 
else
    fprintf('SKIP THIS FILE.\n');return;
end

fprintf('Ignoring selected cases...\n');
ignoreCase = [14;16;23;30];
ignoreSubj = ignoreCase - 6;
if (sum(ismember(caseNum,ignoreCase))>0)
    fprintf('IGNORE THIS CASE: %d\n',caseNum);return;
end

fileName = ['case',num2str(caseNum),'_',cell2mat(routineName(fileIdx)),'_Vol_BR_HR.mat'];
load([dataPath,fileName]);

% Plotting the isovolumetric part in Routine 1
close all
figure('Position',[100,200,600,500])

normIdx = [10, 20].*fs(1);

% In some exceptional case this was wrongly estimated, and creates a
% problem in correlation.

tDevManualNcsThAbd =  0; % Thorax - Abdomen time. 0.53-1.21
if fileNum == 13
    tDevManualNcsThAbd = 0.53-1.21;
end
if fileNum == 43
    tDevManualNcsThAbd = -32.28-(-32.95);
end
if fileNum == 108
    tDevManualNcsThAbd = 3.306-2.9;
end

if abs(tDevManualNcsThAbd) > 0
    fprintf('Correcting time shift...\n');
    nSampDev = uint64(abs(tDevManualNcsThAbd) * fs(1)); % same BIO and NCS frequencies
    % Since we're only looking at isovolumetric part, only making sure ncs
    % thorax and abdomen are correct at the end period, inserting 0s in
    % beginning if needed to correct the time shift, and do not care the slight
    % offset wrt biopac for now.
    if tDevManualNcsThAbd > 0
        ncsRespAbd = [zeros(nSampDev,1);ncsRespAbd(1:end-nSampDev)];
    else
        ncsRespTh = [zeros(nSampDev,1);ncsRespTh(1:end-nSampDev)];
    end
fprintf('Sync: NCS Abd with NCS Th.\n');
end

% NCS sign correction if needed
if fileNum == 33 || fileNum == 108
    ncsRespTh = -ncsRespTh;
end

fSz = 9;
nFig = 2;
plotTst = 320;
ax(1) = subplot(nFig,1,1);
plot(tNcs-plotTst,bio(:,2)./rms(bio(normIdx(1):normIdx(2),2)));
hold on
plot(tNcs-plotTst,bio(:,3)./rms(bio(normIdx(1):normIdx(2),2)));
% xlim([320,375])
legend({'Thorax','Abdomen'},'FontSize',fSz,'Location','southwest');
% xlabel('Time (s)','FontSize',fSz);
ylabel('Normalized BIOPAC','FontSize',fSz)
title(['Isovolumetric Exercise: Case: ',num2str(caseNum)],'FontSize',fSz-2)
% ylim([-4,4])
ax(2) = subplot(nFig,1,2);
plot(tNcs-plotTst,ncsRespTh./rms(ncsRespTh(normIdx(1):normIdx(2))));
hold on
plot(tNcs-plotTst,ncsRespAbd./rms(ncsRespAbd(normIdx(1):normIdx(2))));
% xlim([320,375])
% legend({'Thorax','Abdomen'},'FontSize',fSz,'Location','southwest');
xlabel('Time (s)','FontSize',fSz);ylabel('Normalized NCS','FontSize',fSz)
% title(['Isovolumetric Exercise: Case: ',num2str(caseNum)],'FontSize',fSz-2)
linkaxes(ax,'x');
xlim([0,60]);
% ylim([-4,4])

% Biopac, NCS, Thorax and Abdomen respectively
idxIsoVol = [plotTst+(1/fs(1)),tNcs(end)].*fs;
isoBNTA = [bio(idxIsoVol(1):idxIsoVol(2),2)./rms(bio(normIdx(1):normIdx(2),2)),...  % BIOPAC Thorax
           bio(idxIsoVol(1):idxIsoVol(2),3)./rms(bio(normIdx(1):normIdx(2),3)),...  % BIOPAC Abd
           ncsRespTh(idxIsoVol(1):idxIsoVol(2))./rms(ncsRespTh(normIdx(1):normIdx(2))),... % NCS Thorax
           ncsRespAbd(idxIsoVol(1):idxIsoVol(2))./rms(ncsRespAbd(normIdx(1):normIdx(2)))]; % NCS Abdomen
tIso = tNcs(idxIsoVol(1):idxIsoVol(2))-plotTst;

% =========================================================================
% Keep updating true indicator based on each case. True is marked if at
% least one of NCS and Biopac can see paradoxical movement - so manual
% indicator. 1 if detected, 0 otherwise.

% Expected: [20,23; 35,38; 50,53] Time in sec [start, stop]

[tIsoVolTrue,tCalibIsoVol] = getTimeIsoVolTrue(caseNum);
% tCalibIsoVol = [1,5];
fprintf('Calibration window: [%3.2f, %3.2f].\n',tCalibIsoVol(1),tCalibIsoVol(2));
numIVTrue = size(tIsoVolTrue,1);

%% ------------------------------------------------------------------------
% Developing feature vectors for IsoVolumetric motion detection
% -------------------------------------------------------------------------
slopeBNTA = diff(isoBNTA).*fs(1); 
slopeBNTA = [slopeBNTA(1,:); slopeBNTA]; % Same row size as t

% Rescaling slope, for cases when it shoots up a lot
slopeBNTA = tanh(slopeBNTA);

numFeatIV = 3; 
featDescriptIV = {'Avg SP','StdDev SP','Correlation'};
% Features: [slope product, Avg slope product, Std Dev slope product, correlation] 
% Moving average and standard deviation are taken on a window of length
% tWinFeat, and updated only after tWinSlideFeat
tWinFeat = 1; % Estimate (avg, stddev) features over tWinFeat seconds
tWinCorrFeat = 3; % (s) Window for correlation feat, centered at current index.
tWinSlideFeat = 1; % This is the duration between each window in seconds
numWinFeat = floor((tIso(end)-tWinFeat)/tWinSlideFeat)+1;

slopeProdBio = (slopeBNTA(:,1).*slopeBNTA(:,2));
slopeProdNcs = (slopeBNTA(:,3).*slopeBNTA(:,4));

movAvgBioSP = movavg(slopeProdBio,'Linear',tWinFeat*fs(1)); % Moving Avg over past tWinFeat s
movStdBioSP = movstd(slopeProdBio,[tWinFeat*fs(1)-1,0]); % Moving StdDev over past tWinFeat s
movAvgNcsSP = movavg(slopeProdNcs,'Linear',tWinFeat*fs(1));
movStdNcsSP = movstd(slopeProdNcs,[tWinFeat*fs(1)-1,0]);
corrBioNcsTA = zeros(length(tIso),2);
for iter = 1:length(tIso)
    if iter <= floor(tWinCorrFeat*fs(1)/2)
        winStartIdx = 1;
    else
        winStartIdx = iter - floor(tWinCorrFeat*fs(1)/2);
    end
    winEndIdx = winStartIdx + tWinCorrFeat*fs(1)-1;
    if winEndIdx > length(tIso)
        winEndIdx = length(tIso);
    end
    winCorrBio = corrcoef(isoBNTA(winStartIdx:winEndIdx,1),isoBNTA(winStartIdx:winEndIdx,2));
    corrBioNcsTA(iter,1) = winCorrBio(1,2);
    winCorrNcs = corrcoef(isoBNTA(winStartIdx:winEndIdx,3),isoBNTA(winStartIdx:winEndIdx,4));
    corrBioNcsTA(iter,2) = winCorrNcs(1,2);
end

featBioIV = zeros(numWinFeat,numFeatIV); 
featNcsIV = zeros(numWinFeat,numFeatIV);
tIdxFeatIV = zeros(numWinFeat,1);
for iter = 1:numWinFeat
    % Just keep every ith index of the estimated features, starting from
    % the end of first window index (as moving average and moving standard
    % deviation are estimated with current index and past samples as
    % window.
    idx = (iter)*tWinSlideFeat*fs(1);
    if idx > length(tIso)
        idx = length(tIso);
    end
    featBioIV(iter,1) = movAvgBioSP(idx);
    featBioIV(iter,2) = movStdBioSP(idx);
    featBioIV(iter,3) = corrBioNcsTA(idx,1);

    featNcsIV(iter,1) = movAvgNcsSP(idx);
    featNcsIV(iter,2) = movStdNcsSP(idx);
    featNcsIV(iter,3) = corrBioNcsTA(idx,2);
    tIdxFeatIV(iter) = idx;
end

tIsoFeat = tIso(tIdxFeatIV);
% -------------------------------------------------------------------------
% Plotting features on previous graph
plot(ax(1),tIsoFeat,featBioIV(:,1));
plot(ax(1),tIsoFeat,featBioIV(:,2));
plot(ax(1),tIsoFeat,featBioIV(:,3));
legend(ax(1),{'Thorax','Abdomen','MA','MStdDev','Corr'},'FontSize',fSz,...
    'Location','southwest','Orientation','Horizontal');

plot(ax(2),tIsoFeat,featNcsIV(:,1));
plot(ax(2),tIsoFeat,featNcsIV(:,2));
plot(ax(2),tIsoFeat,featNcsIV(:,3));
% -------------------------------------------------------------------------

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

%% Accuracy estimate
accNcsIV = zeros(numIVTrue,1); % NCS Detected "True" IV corresponding to each True IV
accBioIV = zeros(numIVTrue,1); % BIOPAC Detected "True" IV corresponding to each True IV
% meanNcsMAslopeIV = zeros(numIVTrue,1); % Mean MA slope during each IV
% meanBioMAslopeIV = zeros(numIVTrue,1);

% Only concentrating on windows where IV is detected
for iter = 1:numIVTrue
    indIVTrue = ((tIsoFeat >= tIsoVolTrue(iter,1))&(tIsoFeat <= tIsoVolTrue(iter,2)));
    
%     % Saving mean MA slope during the TRUE IV maneuvers
%     meanNcsMAslopeIV(iter) = mean(featNcsIV(indIVTrue,1));
%     meanBioMAslopeIV(iter) = mean(featBioIV(indIVTrue,1));
%     if meanNcsMAslopeIV(iter)<threshIVNcs(1)
%         accNcsIV(iter) =  1;
%     end
%     if meanBioMAslopeIV(iter)<threshIVBio(1)
%         accBioIV(iter) =  1;
%     end
    if sum(indIVBio(indIVTrue))>0
        accBioIV(iter) = 1;
    end
    if sum(indIVNcs(indIVTrue))>0
        accNcsIV(iter) = 1;
    end

    indIVBioFalse(indIVTrue) = 0;
    indIVNcsFalse(indIVTrue) = 0;
    
end

% Now, for the rest of the duration, what is the mean score : to get
% overall quality indicator of the NCS and BIOPAC signal.
indIVTrue = zeros(numWinFeat,1);
for iter = 1:numIVTrue
    indIVTrue = indIVTrue | ((tIsoFeat >= tIsoVolTrue(iter,1))&(tIsoFeat <= tIsoVolTrue(iter,2)));
end

indNoIV = ~indIVTrue; % Where Iso Vol doesn't exist
tIVBioFalse = sum(indIVBioFalse).*tWinSlideFeat % BIOPAC: Wrong IV estimated time
tIVNcsFalse = sum(indIVNcsFalse).*tWinSlideFeat % NCS: Wrong IV estimated time

% Quality defined as total NON IV time/ time predicted as IV during non IV
% Higher the better, lies in [0,1]
qualIVBio = 1 - (sum(indIVBioFalse)/sum(indNoIV)); 
qualIVNcs = 1 - (sum(indIVNcsFalse)/sum(indNoIV)); 
meanNcsMAslopeNoIV = mean(featNcsIV(indNoIV,1));
meanBioMAslopeNoIV = mean(featBioIV(indNoIV,1));

fprintf('\nCase Num: %d\n',caseNum);
fprintf('True IV events: %d \n',numIVTrue);
fprintf('NCS: IV detected: %d, Quality: %3.2f\n',sum(accNcsIV), qualIVNcs);
fprintf('BIOPAC: IV detected: %d, Quality: %3.2f\n',sum(accBioIV), qualIVBio);

accNcsIV
accBioIV

%%
% close all
figure('Position',[100,50,600,350])

normIdx = [10, 20].*fs(1);

fSz = 8; nFig = 2;
plotTst = 0;
ax(1) = subplot(nFig,1,1);
plot(tIso-plotTst,isoBNTA(:,1));
hold on
plot(tIso-plotTst,isoBNTA(:,2));
plot(tIsoFeat-plotTst,featBioIV(:,1),':','LineWidth',1.5);
plot(tIsoFeat-plotTst,-10.*indIVBio,'-.','LineWidth',1.5);
plot(tIsoFeat-plotTst,-8.*indIVTrue,'-.','LineWidth',1.5);

legend({'Thorax','Abdomen','Slope: MA','Detected Paradoxical Movement'},'FontSize',fSz,'Orientation','Horizontal','EdgeColor',[1,1,1]);
ylabel('Biopac','FontSize',fSz)
% title(['Isovolumetric Exercise: Case: ',num2str(caseNum)],'FontSize',fSz-2)
% ylim([-4,4])
ax(2) = subplot(nFig,1,2);
plot(tIso-plotTst,isoBNTA(:,3));
hold on
plot(tIso-plotTst,isoBNTA(:,4));
plot(tIsoFeat-plotTst,featNcsIV(:,1),':','LineWidth',1.5);
plot(tIsoFeat-plotTst,-10.*indIVNcs,'-.','LineWidth',1.5);
plot(tIsoFeat-plotTst,-8.*indIVTrue,'-.','LineWidth',1.5);

xlabel('Time (s)','FontSize',fSz);
ylabel('NCS','FontSize',fSz)
linkaxes(ax,'x');
xlim([0,60]);
set(ax,'YGrid','on','FontSize',fSz-2);

%%
figure('Position',[100,50,600,350])

normIdx = [10, 20].*fs(1);

fSz = 8; nFig = 2;
plotTst = 0;
ax(1) = subplot(nFig,1,1);
plot(tIso-plotTst,isoBNTA(:,1));
hold on
plot(tIso-plotTst,isoBNTA(:,2)./4);
plot(tIsoFeat-plotTst,4.*indIVBio,':','LineWidth',1.5,'color',[0.47,0.67,0.19]);
plot(tIsoFeat-plotTst,5.*indIVTrue,'--','LineWidth',1);
ylim([-15,10])

legend({'Thorax','Abdomen','Detected','True'},'FontSize',fSz,'Orientation','Horizontal','EdgeColor',[1,1,1]);
ylabel('BIOPAC Respiration','FontSize',fSz)
ax(2) = subplot(nFig,1,2);
plot(tIso-plotTst,isoBNTA(:,3));
hold on
plot(tIso-plotTst,isoBNTA(:,4));
plot(tIsoFeat-plotTst,4.*indIVNcs,':','LineWidth',1.5,'color',[0.47,0.67,0.19]);
plot(tIsoFeat-plotTst,5.*indIVTrue,'--','LineWidth',1);
xlabel('Time (s)','FontSize',fSz);
ylabel('NCS Respiration','FontSize',fSz)
linkaxes(ax,'x');
ylim([-10,10])
xlim([0,60]);
set(ax,'YGrid','on','FontSize',fSz-2);
linkaxes(ax,'x');
xlim([0,50]);
% pause (1)
% close all

%% Save the results
%{
save(['C:\Research\NCS\HumanStudyData\Analysis3_IV3_Feb4\case',num2str(caseNum),'_',cell2mat(routineName(fileIdx)),'_Vol_BR_HR.mat']);

fprintf('Saved the results...\n');
%}


