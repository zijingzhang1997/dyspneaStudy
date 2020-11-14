dataPath=['C:\Sleep test\paper revision\revision\'];
fileName='180SROUTINE1'
fs=10e3;
fsDS=500;   
toff=[10:170]';
%toff=[2:297]';


fsbio=1e3;
fsDSbio=500;
filePathName = [dataPath,fileName,'.tdms'];
convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);
ampTh=ConvertedData.Data.MeasuredData(3).Data;%  rx1 thorax 
bioTh=ConvertedData.Data.MeasuredData(20).Data;  %ch2 biopac thorax 
ampAb=ConvertedData.Data.MeasuredData(8).Data; %rx2 abdomen
bioAb=ConvertedData.Data.MeasuredData(21).Data;  %ch3 biopac abdomen
ampth01=ConvertedData.Data.MeasuredData(13).Data;  %13  tx3
ampAb01=ConvertedData.Data.MeasuredData(18).Data;   %18 tx4

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
ampfiltAb_raw=ampfiltAb;

dsbioTh=resample(bioTh,fsDSbio,fsbio);
dsbioAb=resample(bioAb,fsDSbio,fsbio);
t = ((0:(length(dsbioTh)-1))/fsDSbio)';
dsbioThOff=dsbioTh((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);
dsbioAbOff=dsbioAb((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);

ds_ampth01=resample(ampth01,fsDS,fs);
ds_ampAb01=resample(ampAb01,fsDS,fs);
t = ((0:(length(ds_ampAb01)-1))/fsDSbio)';
ds_ampth01Off=ds_ampth01((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);
ds_ampAb01Off=ds_ampAb01((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);


tOff=((0:(length(dsbioThOff)-1))/fsDSbio)';

bioThfilt = filterLpHp(dsbioThOff,fsDS,opts1); % th amp
bioAbfilt = filterLpHp(dsbioAbOff,fsDS,opts1); 

% opts2.filtType = 'LpHp'; opts2.orderHP = 5;
% opts2.f3db = 0.005; opts2.fpLP = 5; opts2.fstLP = 7;
ds_ampth01filt = filterLpHp(ds_ampth01Off,fsDS,opts1);

ds_ampAb01filt = filterLpHp(ds_ampAb01Off,fsDS,opts1); 
ds_ampAb01filt_raw=ds_ampAb01filt;


[ampfiltAb,PS] = mapminmax(ampfiltAb');
ampfiltAb=ampfiltAb';
[ampfiltTh,PS2] = mapminmax(ampfiltTh');
ampfiltTh=ampfiltTh';
[bioThfilt,PS] = mapminmax(bioThfilt');
bioThfilt=bioThfilt';
[bioAbfilt,PS] = mapminmax(bioAbfilt');
bioAbfilt=bioAbfilt';
[ds_ampth01filt,PS] = mapminmax(ds_ampth01filt');
ds_ampth01filt=ds_ampth01filt';
[ds_ampAb01filt,PS] = mapminmax(ds_ampAb01filt');
ds_ampAb01filt=ds_ampAb01filt';

figure
sz=10;
subplot(6,1,1);
plot(tOff,ampfiltTh,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 220])
title('NCS Th','FontSize',sz)

subplot(6,1,2);
plot(tOff,ampfiltAb,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 168])
title('NCS Ab ','FontSize',sz)

subplot(6,1,3);
plot(tOff,bioThfilt,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 168])
title('chest belt Th','FontSize',sz)

subplot(6,1,4);
plot(tOff,bioAbfilt,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 168])
title('chest belt ab','FontSize',sz)


subplot(6,1,5);
plot(tOff,ds_ampth01filt,'LineWidth',0.5,'color','red');
%plot(ds_ppg,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 220])
title('NCS th notch','FontSize',sz)

subplot(6,1,6);
plot(tOff,ds_ampAb01filt ,'LineWidth',0.5,'color','green');
%plot(ds_ecg,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 220])
title('NCS Ab notch','FontSize',sz)



%%  breath rate estimation  peak-to-peak
opts3.tWinBR = 15; % Window on which br is estimated
opts3.tWin = 2; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.2;

%[brNcs,ncsRespPk] = brEst(ampfiltTh,fsDS,opts3);
[brNcs,ncsRespPk] = brEst(ampfiltAb,fsDS,opts3);
pkMaxNcs = ncsRespPk.idx(ncsRespPk.ind == 1);
pkMinNcs = ncsRespPk.idx(ncsRespPk.ind == 0);
pkMaxNcs = pkMaxNcs(ncsRespPk.idxValidPk);
pkMinNcs = pkMinNcs(ncsRespPk.idxValidPk);

%[brBio,bioRespPk] = brEst(bioThfilt,fsDSbio,opts3);
[brBio,bioRespPk] = brEst(ds_ampAb01filt,fsDS,opts3);
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
xlabel('time (s)')
ylabel('Breath Rate')
legend('chest belt ','NCS ')

%% 


ampfilt1=ampfiltAb_raw;
ampfilt2=ds_ampAb01filt_raw;

[brTh,ThRespPk] = brEst(ampfilt1,fsDS,opts3);
pkMaxTh = ThRespPk.idx(ThRespPk.ind == 1);
pkMinTh = ThRespPk.idx(ThRespPk.ind == 0);
pkMaxTh = pkMaxTh(ThRespPk.idxValidPk);
pkMinTh = pkMinTh(ThRespPk.idxValidPk);
tpkTh=tOff(pkMaxTh);
pkValueTh=ampfilt1(pkMaxTh)-ampfilt1(pkMinTh);
coupleTh=pkValueTh./ampDsOffAb(pkMaxTh);
coupleTh=log2(coupleTh)+11;  %signal quality for wearable NCS 




[brAb,AbRespPk] = brEst(ampfilt2,fsDS,opts3);
pkMaxAb = AbRespPk.idx(AbRespPk.ind == 1);
pkMinAb = AbRespPk.idx(AbRespPk.ind == 0);
pkMaxAb = pkMaxAb(AbRespPk.idxValidPk);
pkMinAb = pkMinAb(AbRespPk.idxValidPk);
tpkAb=tOff(pkMaxAb);
pkValueAb=ampfilt2(pkMaxAb)-ampfilt2(pkMinAb);
coupleAb=pkValueAb./ds_ampAb01Off(pkMaxAb);
coupleAb=log2(coupleAb)+11;  %signal quality for notched NCS 

std_cpTh=std(coupleTh);
mean_cpTh=mean(coupleTh);
std_cpAb=std(coupleAb);
mean_cpAb=mean(coupleAb);

snrTh = snr(ampfilt1,fsDS,3);
snrAb = snr(ampfilt2,fsDS,3);

figure
subplot(4,1,1);
plot(tOff,ampfiltAb,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 168])
title('wearable sensor','FontSize',sz)
subplot(4,1,2);
plot(tOff,-ds_ampAb01filt ,'LineWidth',1,'color','green');
%plot(ds_ecg,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 220])
title('notched sensor','FontSize',sz)



subplot(4,1,3)
plot(tOff,brNcs,'LineWidth',1,'color','red');
hold on
plot(tOff,brBio,'LineWidth',1,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Breath Rate','FontSize',sz)
legend('wearable sensor','notched sensor','FontSize',sz)

subplot(4,1,4)
plot(tpkTh,coupleTh,'color','red','LineWidth',1)
hold on
plot(tpkAb,coupleAb,'color','green','LineWidth',1)
xlabel('time (s)','FontSize',sz)
ylabel('Signal Quality (bits)','FontSize',sz)
