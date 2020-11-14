dataPath=['C:\Sleep test\Data\624\'];
fileName='edwin2'
fs=10e3;
fsDS=500;
toff=[10:230]';

fsbio=1e3;
fsDSbio=500;
filePathName = [dataPath,fileName,'.tdms'];
convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);
ampTh=ConvertedData.Data.MeasuredData(8).Data;
bioTh=ConvertedData.Data.MeasuredData(20).Data;
ampAb=ConvertedData.Data.MeasuredData(3).Data;
bioAb=ConvertedData.Data.MeasuredData(21).Data;


ampdsTh=resample(ampTh,fsDS,fs);
ampdsAb=resample(ampAb,fsDS,fs);
t = ((0:(length(ampdsTh)-1))/fsDS)';
ampDsOffTh=ampdsTh((toff(1)*fsDS):toff(size(toff))*fsDS);
ampDsOffAb=ampdsAb((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOffTh)-1))/fsDS)';
opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 1; opts1.fstLP = 1.5;
ampfiltTh = filterLpHp(ampDsOffTh,fsDS,opts1); % th amp
ampfiltAb = filterLpHp(ampDsOffAb,fsDS,opts1); 

dsbioTh=resample(bioTh,fsDSbio,fsbio);
dsbioAb=resample(bioAb,fsDSbio,fsbio);
t = ((0:(length(dsbioTh)-1))/fsDSbio)';
dsbioThOff=dsbioTh((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);
dsbioAbOff=dsbioAb((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);

tOff=((0:(length(dsbioThOff)-1))/fsDSbio)';

% opts1.filtType = 'LpHp'; opts1.orderHP = 5;
% opts1.f3db = 0.05; opts1.fpLP = 5; opts1.fstLP = 7;
bioThfilt = filterLpHp(dsbioThOff,fsDS,opts1); % th amp
bioAbfilt = filterLpHp(dsbioAbOff,fsDS,opts1); 


opts2.filtType = 'LpHp'; opts2.orderHP = 5;
opts2.f3db = 0.7; opts2.fpLP = 1.5; opts2.fstLP = 2;
ampfiltHRTh = filterLpHp(ampDsOffTh,fsDS,opts2); 

[ampfiltAb,PS] = mapminmax(ampfiltAb');
ampfiltAb=ampfiltAb';
[ampfiltTh,PS2] = mapminmax(ampfiltTh');
ampfiltTh=ampfiltTh';
[bioThfilt,PS] = mapminmax(bioThfilt');
bioThfilt=bioThfilt';
[bioAbfilt,PS] = mapminmax(bioAbfilt');
bioAbfilt=bioAbfilt';
figure
sz=10;
subplot(4,1,1);
plot(tOff,ampfiltTh,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])
title('NCS Th 2.4GHz','FontSize',sz)

subplot(4,1,2);
plot(tOff,ampfiltAb,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])
title('NCS Ab 2.4GHz','FontSize',sz)

subplot(4,1,3);
plot(tOff,bioThfilt,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])
 title('BIOPAC Th 2.4GHz','FontSize',sz)

subplot(4,1,4);
plot(tOff,bioAbfilt,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])
title('BIOPAC Ab 2.4GHz','FontSize',sz)





%%  breath rate estimation  peak-to-peak
opts3.tWinBR = 15; % Window on which br is estimated
opts3.tWin = 4; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.3;

[brNcs,ncsRespPk] = brEst(ampfiltAb,fsDS,opts3);
pkMaxNcs = ncsRespPk.idx(ncsRespPk.ind == 1);
pkMinNcs = ncsRespPk.idx(ncsRespPk.ind == 0);
pkMaxNcs = pkMaxNcs(ncsRespPk.idxValidPk);
pkMinNcs = pkMinNcs(ncsRespPk.idxValidPk);

[brBio,bioRespPk] = brEst(bioAbfilt,fsDSbio,opts3);
pkMaxBio = bioRespPk.idx(bioRespPk.ind == 1);
pkMinBio = bioRespPk.idx(bioRespPk.ind == 0);
pkMaxBio = pkMaxBio(bioRespPk.idxValidPk);
pkMinBio = pkMinBio(bioRespPk.idxValidPk);

BRdiff=brBio-brNcs;
% rmseBR = sqrt(mean((brBio-brNcs).^2));
% fprintf('Done. RMSE BR: %3.2f\n',rmseBR);

figure()
plot(tOff,brBio,'LineWidth',0.5,'color','red');

hold on
plot(tOff,brNcs,'LineWidth',0.5,'color','green');
xlabel('t/s')
ylabel('Breath Rate')
legend('Biopac Ab','amp NCS Ab')
%% bA plot correlation plot
pkLS = regress(brBio,[ones(length(brNcs),1) brNcs]);
pkBaMean = mean(brBio - brNcs);
pkBaStd = std(brBio - brNcs);
pkBaStdLim = [pkBaStd*1.96+pkBaMean, -pkBaStd*1.96+pkBaMean]; 
[rPkVol,pPkVol] = corrcoef(brNcs,brBio);
mean=(brNcs+brBio)./2;
diff=brNcs-brBio;

a=[4,61,76,162,195,200];sz=8;% for 0.9 side R
figure()
scatter(brNcs(a(1)*fsDS:a(2)*fsDS),brBio(a(1)*fsDS:a(2)*fsDS),'o','filled'); % make into different parts
hold on
scatter(brNcs(a(2)*fsDS:a(3)*fsDS),brBio(a(2)*fsDS:a(3)*fsDS),'+');
scatter(brNcs(a(3)*fsDS:a(4)*fsDS),brBio(a(3)*fsDS:a(4)*fsDS),'*');
scatter(brNcs(a(4)*fsDS:a(5)*fsDS),brBio(a(4)*fsDS:a(5)*fsDS),'x');
scatter(brNcs(a(5)*fsDS:a(6)*fsDS),brBio(a(5)*fsDS:a(6)*fsDS),'o','filled');
legend('normal breath','stop breath','deep breath','fast breath','normal breath 2');
xVol = 0:10:60;
plot(xVol,pkLS(1)+pkLS(2)*xVol,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
text(50,30,['\bf', 'r=',num2str(rPkVol(1,2),3)],'FontSize',11)
xlabel('NCS breath rate(BPM)','FontName',font,'FontSize',sz);

figure()
scatter(mean(a(1)*fsDS:a(2)*fsDS),diff(a(1)*fsDS:a(2)*fsDS),'o'); % Normal
xVol = 0:10:60;
hold on
scatter(mean(a(2)*fsDS:a(3)*fsDS),diff(a(2)*fsDS:a(3)*fsDS),'+');
scatter(mean(a(3)*fsDS:a(4)*fsDS),diff(a(3)*fsDS:a(4)*fsDS),'*');
scatter(mean(a(4)*fsDS:a(5)*fsDS),diff(a(4)*fsDS:a(5)*fsDS),'x');
scatter(mean(a(5)*fsDS:a(6)*fsDS),diff(a(5)*fsDS:a(6)*fsDS),'o','filled');
plot(xVol,pkBaMean.*ones(length(xVol),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
plot(xVol,pkBaStdLim(1).*ones(length(xVol),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
plot(xVol,pkBaStdLim(2).*ones(length(xVol),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
ylim([-40 40]);xlim([0 60]);
legend('normal breath','stop breath','deep breath','fast breath','normal breath 2');

text(45,pkBaMean+5,['\bf', num2str(pkBaMean,3)],'FontSize',11)
text(45,pkBaStdLim(1)+5,['\bf', num2str(pkBaStdLim(1),3)],'FontSize',11)
text(45,pkBaStdLim(2)+5,['\bf', num2str(pkBaStdLim(2),3)],'FontSize',11)
%%
% heart beat rate filter
% yampDsOff=fft(ampfiltHRTh);
% f = (0:length(yampDsOff)-1)*fsDS/length(yampDsOff)';
% yam=abs(yampDsOff);
% f=f(1,1:1000);
% yam=yam(1:1000,1)';
% yammax=max(yam);
% yammin=min(yam);
% yamnorm=zeros(1,1000);
% yamnorm=(yam-yammin)/(yammax-yammin);
% figure
% plot(f,yamnorm,'LineWidth',0.5,'color','blue');
% xlabel('frequency/Hz')
% 
% title('amp Th')


% % %%  Using cross correlation to estimate time shift, in case
% tCorr = [5, 50]; 
% fprintf('\nTime sync, Calibration: correlating window [%d,%d]s\n',tCorr(1),tCorr(2));
% 
% nStart = tCorr(1)*fsDS; % Assuming same fs for both ncs and biopac
% nEnd = tCorr(2)*fsDS; 
% maxLag = 1000;
% [r, lags] = xcorr(ampfiltThnorm(nStart:nEnd),ampfiltAbnorm(nStart:nEnd),maxLag);
% figure; plot(lags,r)
% [~,rMaxIdx] = max(abs(r));
% lagsMax = lags(rMaxIdx);
% tDevCalib = lagsMax/fsDS;
% fprintf('Suggested NCS Th-Bio Th calibration time offset is %f\n',tDevCalib);   
%  
% 
% 
% %% ONLY call this block once!!
% % Now shift the waveforms to compensate the time difference
% % Only shift if time difference is more than a threshold and less than a
% % threshold: Ideally minimum should be sampling frequency
% 
% fprintf('Performing synchronization based on time-shift estimate...  \n ');
% thMinCorr = 0.4;
% tOffMinMax = [0.005, 1.0];
% % First the calibration waveforms (Bio, NCS Th):
% if abs(tDevCalib)>=tOffMinMax(1) && abs(tDevCalib) <= tOffMinMax(2) 
%     nSampDev = abs(tDevCalib) * fsDS; % same BIO and NCS frequencies
%     if tDevCalib > 0
%         ampfiltThnorm = ampfiltThnorm(nSampDev+1:end,:);
%         ampfiltAbnorm = ampfiltAbnorm(1:end-nSampDev,:);
%     else
%         ampfiltThnorm = ampfiltThnorm(1:end-nSampDev,:);
%         ampfiltAbnorm = ampfiltAbnorm(nSampDev+1:end,:);
%     end
%     toffcorr = tOff(1:end-nSampDev);
%     fprintf('\nSync: NCS calib Th with Biopac.\n');
% else
%     tDevCalib = 0; % To not be recorded, since we are not correcting
% end
% 
% figure()
% plot(toffcorr,ampfiltThnorm,'LineWidth',0.5,'color','red');
% 
% hold on
% plot(toffcorr,ampfiltAbnorm,'LineWidth',0.5,'color','green');
% xlabel('t/s')
% ylabel('normalized Amplitude')
% legend('amp NCS TH','amp NCS AB')