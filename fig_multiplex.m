dataPath=['F:\Sleep test\Data\414\'];
fileName='0.9g1.9g';
fs=10e3;
fsDS=500;
toff=[10:170]';
sz=11;
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
[ampfilt1,PS] = mapminmax(ampfilt1');
ampfilt1=ampfilt1';
[ampfilt2,PS] = mapminmax(ampfilt2');
ampfilt2=ampfilt2';

figure()

subplot(2,1,1);
plot(tOff,ampfilt1,'LineWidth',0.5,'color','red');
xlabel('time(s)','FontSize',sz)
ylabel('Amp(a.u.)','FontSize',sz)

 title('0.9Ghz multiplexing','FontSize',sz)
 
subplot(2,1,2);
plot(tOff,ampfilt2,'LineWidth',0.5,'color','blue');
xlabel('time(s)','FontSize',sz)
ylabel('Amp(a.u.)','FontSize',sz)

title('1.5Ghz multiplexing','FontSize',sz)


 
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
coupleTh=100*coupleTh;




[brAb,AbRespPk] = brEst(ampfilt2,fsDS,opts3);
pkMaxAb = AbRespPk.idx(AbRespPk.ind == 1);
pkMinAb = AbRespPk.idx(AbRespPk.ind == 0);
pkMaxAb = pkMaxAb(AbRespPk.idxValidPk);
pkMinAb = pkMinAb(AbRespPk.idxValidPk);
tpkAb=tOff(pkMaxAb);
pkValueAb=ampfilt2(pkMaxAb)-ampfilt2(pkMinAb);
coupleAb=pkValueAb./ampDsOff2(pkMaxAb);
coupleAb=100*coupleAb;
%  figure()
%  plot(tpkAb,pkValueAb)
figure()

subplot(3,1,1);
plot(tOff,ampfilt1,'LineWidth',0.5,'color','red');
xlabel('time(s)','FontSize',sz)
ylabel('Amp(a.u.)','FontSize',sz)


 
subplot(3,1,2);
plot(tOff,ampfilt2,'LineWidth',0.5,'color','blue');
xlabel('time(s)','FontSize',sz)
ylabel('Amp(a.u.)','FontSize',sz)



subplot(3,1,3);
plot(tOff,brTh,'LineWidth',1,'color','red');

hold on
plot(tOff,brAb,'LineWidth',1,'color','blue');
xlabel('time (s)','FontSize',sz)
legend('0.9GHz','1.5GHz','FontSize',sz)
ylabel('BR (BPM)','FontSize',sz) 


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
ampfilt1norm=ampfilt1norm+mean(ampfilt2norm)-mean(ampfilt1norm);

% figure()
% dtw(ampfilt1norm,ampfilt2norm);
% [dist,ix,iy] = dtw(ampfilt1norm,ampfilt2norm);
% figure()
% plot(tOff,ampfilt1norm,'LineWidth',0.5,'color','red');
% 
% hold on
% plot(tOff,ampfilt2norm,'LineWidth',0.5,'color','green');
% xlabel('t/s')
% ylabel('normalized Amplitude')
% legend('amp 2.4GHz','amp 1.8GHz')


% figure()
% plot(ampfilt1norm(ix))
% hold on
% plot(ampfilt2norm(iy))

% figure()
% plot(tpkTh,coupleTh)
% hold on
% plot(tpkAb,coupleAb)
% xlabel('t/s')
% ylabel('coupling strength(%)')
% legend('th','ab')
