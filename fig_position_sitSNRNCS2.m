dataPath=['F:\Sleep test\Data\519\'];
fileName='50k'
fs=10e3;
fsDS=100;
toff=[30:120]';



filePathName = [dataPath,fileName,'.tdms'];
convertTDMS(true,filePathName);


load([dataPath,fileName,'.mat']);
ampTh=ConvertedData.Data.MeasuredData(3).Data;

ampAb=ConvertedData.Data.MeasuredData(4).Data;


ampdsTh=resample(ampTh,fsDS,fs);
ampdsAb=resample(ampAb,fsDS,fs);
t = ((0:(length(ampdsTh)-1))/fsDS)';
ampfiltTh=ampdsTh((toff(1)*fsDS):toff(size(toff))*fsDS);
ampfiltAb=ampdsAb((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampfiltTh)-1))/fsDS)';
opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 1; opts1.fstLP = 1.5;
ampfiltTh = filterLpHp(ampfiltTh,fsDS,opts1); % th amp
ampfiltAb = filterLpHp(ampfiltAb,fsDS,opts1); 



tOff=((0:(length(ampfiltTh)-1))/fsDS)';

% opts1.filtType = 'LpHp'; opts1.orderHP = 5;
% opts1.f3db = 0.05; opts1.fpLP = 5; opts1.fstLP = 7;





[ampfiltAb,PS] = mapminmax(ampfiltAb');
ampfiltAb=ampfiltAb';
[ampfiltTh,PS] = mapminmax(ampfiltTh');
ampfiltTh=ampfiltTh';


figure
sz=10;
subplot(2,1,1);
plot(tOff,ampfiltTh,'LineWidth',0.5,'color','red');
xlabel('time(s)','FontSize',sz)
ylabel('Amp(a.u.)','FontSize',sz)
xlim([0 50])

title('Th recumbent','FontSize',sz)

subplot(2,1,2);
plot(tOff,ampfiltAb,'LineWidth',0.5,'color','green');
xlabel('time(s)','FontSize',sz)
ylabel('Amp(a.u.)','FontSize',sz)
xlim([0 50])


title('Abd recumbent','FontSize',sz)

SNRampfiltTh = snr(ampfiltTh,fsDS,2);
SNRampfiltAb = snr(ampfiltAb,fsDS,2);


figure()
snr(ampfiltAb,fsDS,3);



[r,noisepow] =snr(ampfiltAb,fsDS,3);