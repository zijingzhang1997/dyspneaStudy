clear all
dataPath=['C:\Sleep test\Data\624\'];
fileName='pragya';
fs=10e3;
fsDS=500;
toff=[10:230]';
sz=10;
toff2=[110,130];  % BR calculation delete isovolumetric [120,140]
fsbio=1e3;
fsDSbio=500;
filePathName = [dataPath,fileName,'.tdms'];
%convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);


ampTh=ConvertedData.Data.MeasuredData(3).Data; %tx1rx1
bioTh=ConvertedData.Data.MeasuredData(20).Data;
ampAb=ConvertedData.Data.MeasuredData(8).Data;  %tx2rx2
bioAb=ConvertedData.Data.MeasuredData(21).Data;


ampdsTh=resample(ampTh,fsDS,fs);
ampdsAb=resample(ampAb,fsDS,fs);
t = ((0:(length(ampdsTh)-1))/fsDS)';
ampDsOffTh=ampdsTh((toff(1)*fsDS):toff(size(toff))*fsDS);
ampDsOffAb=ampdsAb((toff(1)*fsDS):toff(size(toff))*fsDS);
ampDsOffTh=cat(1,ampDsOffTh(1:toff2(1)*fsDS),ampDsOffTh(toff2(2)*fsDS:length(ampDsOffTh)));
ampDsOffAb=cat(1,ampDsOffAb(1:toff2(1)*fsDS),ampDsOffAb(toff2(2)*fsDS:length(ampDsOffAb)));
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
dsbioThOff=cat(1,dsbioThOff(1:toff2(1)*fsDS),dsbioThOff(toff2(2)*fsDS:length(dsbioThOff)));
dsbioAbOff=cat(1,dsbioAbOff(1:toff2(1)*fsDS),dsbioAbOff(toff2(2)*fsDS:length(dsbioAbOff)));
tOff=((0:(length(dsbioThOff)-1))/fsDSbio)';

% opts1.filtType = 'LpHp'; opts1.orderHP = 5;
% opts1.f3db = 0.05; opts1.fpLP = 5; opts1.fstLP = 7;
bioThfilt = filterLpHp(dsbioThOff,fsDS,opts1); % th amp
bioAbfilt = filterLpHp(dsbioAbOff,fsDS,opts1); 


opts2.filtType = 'LpHp'; opts2.orderHP = 5;
opts2.f3db = 0.7; opts2.fpLP = 1.5; opts2.fstLP = 2;
ampfiltHRTh = filterLpHp(ampDsOffTh,fsDS,opts2); 



figure

subplot(4,1,1);
plot(tOff,ampfiltTh,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)

 title('NCS Th ','FontSize',sz)

subplot(4,1,2);
plot(tOff,ampfiltAb,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
title('NCS Ab ')

subplot(4,1,3);
plot(tOff,bioThfilt,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
 title('BIOPAC Th ')

subplot(4,1,4);
plot(tOff,bioAbfilt,'LineWidth',0.5,'color','green');
xlabel('t/s')
title('BIOPAC Ab ')





%%  breath rate estimation  peak-to-peak
opts3.tWinBR = 15; % Window on which br is estimated
opts3.tWin = 2; % Window for peak detection moving average
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

zero_long=zeros(toff2(2)*fsDS-toff2(1)*fsDS,1);
brBio_long=cat(1,brBio(1:toff2(1)*fsDS),zero_long,brBio(toff2(1)*fsDS +1:length(brBio)));
brNcs_long=cat(1,brNcs(1:toff2(1)*fsDS),zero_long,brNcs(toff2(1)*fsDS +1:length(brNcs)));
tOfflong=((0:(length(brNcs_long)-1))/fsDS)';
figure()
plot(tOfflong,brBio_long,'LineWidth',1,'color','red');

hold on
plot(tOfflong,brNcs_long,'LineWidth',1,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('BR (BPM)','FontSize',sz)
xlim([0 220])
legend('BIOPAC Ab','NCS Ab','FontSize',sz)

%% bA plot correlation plot  brBio_1 is without outliers (for calculate mean and std), brBio after making outlier points to zero(for plot)
pkBaMean = mean(brBio - brNcs);
pkBaStd = std(brBio - brNcs);
brBio_1=[];
brNcs_1=[];
for i=1:length(brBio)
    diff=brNcs(i)-brBio(i);
   if diff > 3*pkBaStd || diff< -3*pkBaStd
       brBio(i)=0;brNcs(i)=0;
   else
        brBio_1(end+1)=brBio(i);brNcs_1(end+1)=brNcs(i);
       
   end
end
brBio_1=brBio_1';
brNcs_1=brNcs_1';
ratio=(length(brBio)-length(brBio_1))/length(brBio);

pkLS = regress(brBio_1,[ones(length(brNcs_1),1) brNcs_1]);
pkBaMean = mean(brBio_1 - brNcs_1)
pkBaStd = std(brBio_1 - brNcs_1)
pkBaStdLim = [pkBaStd*1.96+pkBaMean, -pkBaStd*1.96+pkBaMean]; 
[rPkVol,pPkVol] = corrcoef(brNcs_1,brBio_1);
mean=(brNcs+brBio)./2;
diff=brNcs-brBio;

a=[5.8,62.7,75,164,194.5,length(brNcs)/fsDS];% for 1.8
%a=[4,61,76,162,195,200];% for 0.9 side R
figure
subplot(1,2,1);
scatter(brNcs(a(1)*fsDS:a(2)*fsDS),brBio(a(1)*fsDS:a(2)*fsDS),'o','filled'); % make into different parts
hold on
scatter(brNcs(a(2)*fsDS:a(3)*fsDS),brBio(a(2)*fsDS:a(3)*fsDS),'+');
scatter(brNcs(a(3)*fsDS:a(4)*fsDS),brBio(a(3)*fsDS:a(4)*fsDS),'*');
scatter(brNcs(a(4)*fsDS:a(5)*fsDS),brBio(a(4)*fsDS:a(5)*fsDS),'x');
scatter(brNcs(a(5)*fsDS:a(6)*fsDS),brBio(a(5)*fsDS:a(6)*fsDS),'o','filled');
%legend('normal breath','stop breath','deep breath','fast breath','normal breath 2');
xVol = 0:10:60;
plot(xVol,pkLS(1)+pkLS(2)*xVol,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
text(50,30,['\bf', 'r=',num2str(rPkVol(1,2),3)],'FontSize',11)
xlabel('NCS BR (BPM)','FontSize',sz);
ylabel('BIOPAC BR (BPM)','FontSize',sz);
ylim([0 60]);

subplot(1,2,2);
scatter(mean(a(1)*fsDS:a(2)*fsDS),diff(a(1)*fsDS:a(2)*fsDS),'o','filled'); % Normal
xVol = 0:10:60;
hold on
scatter(mean(a(2)*fsDS:a(3)*fsDS),diff(a(2)*fsDS:a(3)*fsDS),'+');
scatter(mean(a(3)*fsDS:a(4)*fsDS),diff(a(3)*fsDS:a(4)*fsDS),'*');
scatter(mean(a(4)*fsDS:a(5)*fsDS),diff(a(4)*fsDS:a(5)*fsDS),'x');
scatter(mean(a(5)*fsDS:a(6)*fsDS),diff(a(5)*fsDS:a(6)*fsDS),'o','filled');
plot(xVol,pkBaMean.*ones(length(xVol),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
plot(xVol,pkBaStdLim(1).*ones(length(xVol),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
plot(xVol,pkBaStdLim(2).*ones(length(xVol),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
xlabel('mean BR (BPM)','FontSize',sz);
ylabel('difference BR (BPM)','FontSize',sz);
ylim([-40 40]);xlim([0 60]);
legend('normal breath','hold breath','deep breath','fast breath','normal breath 2','NumColumns',3);

text(45,pkBaMean+5,['\bf', num2str(pkBaMean,3)],'FontSize',10)
text(45,pkBaStdLim(1)+5,['\bf', num2str(pkBaStdLim(1),3)],'FontSize',10)
text(45,pkBaStdLim(2)+5,['\bf', num2str(pkBaStdLim(2),3)],'FontSize',10)

%% 

%save('C:\Sleep test\Data\624\final\pragya','brBio_1','brNcs_1','brBio','brNcs');
