dataPath=['C:\Sleep test\Data\624\'];
fileName='thomas_recumbent';
fs=10e3;
% fs=1e3; bippac
fsDS=500;
toff=[20:200]';
filePathName = [dataPath,fileName,'.tdms'];
convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);
amp1=ConvertedData.Data.MeasuredData(3).Data; %biopac

amp2=ConvertedData.Data.MeasuredData(8).Data; %biopac




ampds1=resample(amp1,fsDS,fs);

t = ((0:(length(ampds1)-1))/fsDS)';
ampDsOff1=ampds1((toff(1)*fsDS):toff(size(toff))*fsDS);

tOff=((0:(length(ampDsOff1)-1))/fsDS)';
opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 0.8; opts1.fstLP = 1;
ampfilt1 = filterLpHp(ampDsOff1,fsDS,opts1); % th amp


ampds2=resample(amp2,fsDS,fs);

t = ((0:(length(ampds2)-1))/fsDS)';
ampDsOff2=ampds2((toff(1)*fsDS):toff(size(toff))*fsDS);

tOff=((0:(length(ampDsOff2)-1))/fsDS)';

opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 0.8 ;opts1.fstLP = 1;
ampfilt2 = filterLpHp(ampDsOff2,fsDS,opts1); % th amp
%ampfilt2=-ampfilt2;

% figure
% 
% subplot(4,1,1);
% plot(tOff,ampfilt1,'LineWidth',0.5,'color','red');
% xlabel('t/s')
% 
%  title('Amp TH')
% 
% subplot(4,1,2);
% plot(tOff,phfilt1,'LineWidth',0.5,'color','green');
% xlabel('t/s')
% title('Ph TH')
% 
% subplot(4,1,3);
% plot(tOff,ampfilt2,'LineWidth',0.5,'color','red');
% xlabel('t/s')
%  title('Amp AB')
% 
% subplot(4,1,4);
% plot(tOff,phfilt2,'LineWidth',0.5,'color','green');
% xlabel('t/s')
% title('Ph AB')

ampfilt1max=max(ampfilt1);
ampfilt1min=min(ampfilt1);
ampfilt1norm=zeros(1,length(ampfilt1));
ampfilt1norm=(ampfilt1-ampfilt1min)/(ampfilt1max-ampfilt1min);


ampfilt2max=max(ampfilt2);
ampfilt2min=min(ampfilt2);
ampfilt2norm=zeros(1,length(ampfilt2));
ampfilt2norm=(ampfilt2-ampfilt2min)/(ampfilt2max-ampfilt2min);



%xlim([157 217])


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
tWinFeat = 5; % Estimate (avg, stddev) features over tWinFeat seconds
tWinCorrFeat = 5; % (s) Window for correlation feat, centered at current index.
tWinSlideFeat = 5; % This is the duration between each window in seconds
numWinFeat = floor((tOff(end)-tWinFeat)/tWinSlideFeat)+1;


slopeProdNcs = (slopeBNTA1.*slopeBNTA2)

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
t_iso=linspace(0,270,2700);
iso=zeros(length(t_iso),1);
% iso_ind=[80,100,120,140,160,180,200,220];
iso_ind=[110,130];
iso_ind=iso_ind*10;
iso(iso_ind(1):iso_ind(2))=-1;

slope=featNcsIV(:,1);
slope(slope<0)=-1;
slope(slope>0)=0;

figure()
sz=10;
subplot(2,1,1);
plot(tOff,ampfilt1norm,'LineWidth',0.8,'color','red');

hold on
plot(tOff,ampfilt2norm,'LineWidth',0.8,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
legend('Th','Ab','FontSize',sz)
 xlim([20 200])
subplot(2,1,2);
plot(tIsoFeat,slope,'LineWidth',2);
hold on 
plot(t_iso,iso,':','LineWidth',3);
xlabel('time (s)','FontSize',sz)
title('isovolumetric detection','FontSize',sz)
xlim([20 200])
legend('detected','true','FontSize',sz)

% -------------------------------------------------------------------------
