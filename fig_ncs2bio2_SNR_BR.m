dataPath=['F:\Sleep test\Data\317NCS2bio2\'];
fileName='1.8R'
fs=10e3;
fsDS=500;
toff=[10:55]';

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

SNRampfiltTh = snr(ampfiltTh,fsDS);
SNRampfiltAb = snr(ampfiltAb,fsDS);
SNRbioThfilt= snr(bioThfilt,fsDS);
SNRbioAbfilt= snr(bioAbfilt,fsDS);

subplot(4,1,1);
plot(tOff,ampfiltTh,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])
%title('NCS Th supine','FontSize',sz)

subplot(4,1,2);
plot(tOff,ampfiltAb,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])
%title('NCS Ab supine','FontSize',sz)

subplot(4,1,3);
plot(tOff,bioThfilt,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])
 %title('BIOPAC Th supine','FontSize',sz)

subplot(4,1,4);
plot(tOff,bioAbfilt,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])
%title('BIOPAC Ab supine','FontSize',sz)




