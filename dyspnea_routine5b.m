close all
dataPath=['C:\Sleep test\dyspnea\data\10_6\'];
fileName='Case1Routine5b';
fs=10e3;
fsDS=500;
toff=[1:178]';
period=[toff(1),toff(end)];
tNoise=[];
%toff=[10:280]';

fsbio=1e3;
fsDSbio=500;
filePathName = [dataPath,fileName,'.tdms'];
filePathName_m = [dataPath,fileName,'.mat'];
if ~exist(filePathName_m,'file')
   convertTDMS(true,filePathName);
end
load(filePathName_m);


load([dataPath,fileName,'.mat']);
ampTh=ConvertedData.Data.MeasuredData(3).Data;%  rx1 thorax 
bioTh=ConvertedData.Data.MeasuredData(20).Data;  %ch2 biopac thorax 
ampTh_no=ConvertedData.Data.MeasuredData(13).Data;  %13  tx3
ampAb_no=ConvertedData.Data.MeasuredData(18).Data;   %18 tx4
ampAb=ConvertedData.Data.MeasuredData(8).Data; %rx2 abdomen
bioAb=ConvertedData.Data.MeasuredData(21).Data;  %ch3 biopac abdomen
PPG=ConvertedData.Data.MeasuredData(22).Data;  %19 biopac ch1
ECG=ConvertedData.Data.MeasuredData(19).Data;



ampdsTh=resample(ampTh,fsDS,fs);
ampdsAb=resample(ampAb,fsDS,fs);

ampDsOffTh=ampdsTh((toff(1)*fsDS):toff(size(toff))*fsDS);
ampDsOffAb=ampdsAb((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampDsOffTh)-1))/fsDS)';
opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 1; opts1.fstLP = 1.2;
ampfiltTh = filterLpHp(ampDsOffTh,fsDS,opts1); % th amp
ampfiltAb = filterLpHp(ampDsOffAb,fsDS,opts1); 

ampdsTh_no=resample(ampTh_no,fsDS,fs);
ampdsAb_no=resample(ampAb_no,fsDS,fs);

ampDsOffTh_no=ampdsTh_no((toff(1)*fsDS):toff(size(toff))*fsDS);
ampDsOffAb_no=ampdsAb_no((toff(1)*fsDS):toff(size(toff))*fsDS);
ampfiltTh_no = filterLpHp(ampDsOffTh_no,fsDS,opts1); % th amp
ampfiltAb_no = filterLpHp(ampDsOffAb_no,fsDS,opts1); 


dsbioTh=resample(bioTh,fsDSbio,fsbio);
dsbioAb=resample(bioAb,fsDSbio,fsbio);

dsbioThOff=dsbioTh((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);
dsbioAbOff=dsbioAb((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);

ds_ppg=resample(PPG,fsDSbio,fsbio);
ds_ecg=resample(ECG,fsDSbio,fsbio);

ds_ppgOff=ds_ppg((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);
ds_ecgOff=ds_ecg((toff(1)*fsDSbio):toff(size(toff))*fsDSbio);


tOff=((0:(length(dsbioThOff)-1))/fsDSbio)';

bioThfilt = filterLpHp(dsbioThOff,fsDS,opts1); % th amp
bioAbfilt = filterLpHp(dsbioAbOff,fsDS,opts1); 

opts5.filtType = 'LpHp';
opts5.f3db = 4; opts5.fpLP = 20; opts5.fstLP = 25;
% opts2.filtType = 'LpHp'; opts2.orderHP = 5;
% opts2.f3db = 0.005; opts2.fpLP = 5; opts2.fstLP = 7;
% ds_ppgfilt = filterLpHp(ds_ppgOff,fsDS,opts2); 
ds_ecgfilt = filterLpHp(ds_ecgOff,fsDS,opts5); 
%ds_ecgfilt=ds_ecgOff;
ds_ppgfilt=ds_ppgOff;


[ampfiltAb_norm,PS] = mapminmax(ampfiltAb');
ampfiltAb_norm=ampfiltAb_norm';
[ampfiltTh_norm,PS2] = mapminmax(ampfiltTh');
ampfiltTh_norm=ampfiltTh_norm';
[ampfiltAb_no_norm,PS] = mapminmax(ampfiltAb_no');
ampfiltAb_no_norm=ampfiltAb_no_norm';
[ampfiltTh_no_norm,PS2] = mapminmax(ampfiltTh_no');
ampfiltTh_no_norm=ampfiltTh_no_norm';
[bioThfilt_norm,PS] = mapminmax(bioThfilt');
bioThfilt_norm=bioThfilt_norm';
[bioAbfilt_norm,PS] = mapminmax(bioAbfilt');
bioAbfilt_norm=bioAbfilt_norm';
[ds_ppgfilt_norm,PS] = mapminmax(ds_ppgfilt');
ds_ppgfilt_norm=ds_ppgfilt_norm';
[ds_ecgfilt_norm,PS] = mapminmax(ds_ecgfilt');
ds_ecgfilt_norm=ds_ecgfilt_norm';

h(1)=figure;
sz=10;
subplot(6,1,1);
plot(tOff,ampfiltTh_norm,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 220])

titleName=[fileName,'  NCS Th'];
title(titleName,'FontSize',sz);

subplot(6,1,2);
plot(tOff,ampfiltAb_norm,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 168])
title('NCS Ab ','FontSize',sz)

subplot(6,1,3);
plot(tOff,ampfiltTh_no_norm,'LineWidth',0.5,'color','red');
%plot(ds_ppg,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 220])
title('notch th','FontSize',sz)

subplot(6,1,4);
plot(tOff,ampfiltAb_no_norm,'LineWidth',0.5,'color','green');
%plot(ds_ecg,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 220])
title('notch Ab','FontSize',sz)

subplot(6,1,5);
plot(tOff,bioThfilt_norm,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 168])
title('chest belt Th','FontSize',sz)

subplot(6,1,6);
plot(tOff,bioAbfilt_norm,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 168])
title('chest belt ab','FontSize',sz)

set(gcf,'Position',[100,200,1000,800]);
figName = [dataPath,'fig\',fileName,'waveform'];


% figure()
% subplot(2,1,1);
% plot(tOff,ds_ppgfilt,'LineWidth',0.5,'color','red');
% xlabel('time (s)','FontSize',sz)
% ylabel('Amp (a.u.)','FontSize',sz)
% xlim([0 20])
% titleName=[fileName,'  PPG']
% title(titleName,'FontSize',sz)
% 
% 
% subplot(2,1,2);
% plot(tOff,ds_ecgfilt,'LineWidth',0.5,'color','green');
% xlabel('time (s)','FontSize',sz)
% ylabel('Amp (a.u.)','FontSize',sz)
% xlim([0 20])
% title('ECG','FontSize',sz)
% 
% 
% set(gcf,'Position',[100,200,500,350]);
% figName = [dataPath,'fig/',fileName,'Ecg_ppg'];
% %  saveas(gca,figName,'tif')
% % saveas(gca,figName,'fig')

%%  breath rate estimation  peak-to-peak

%ampfilt1=ampfiltTh_norm;
%ampfilt2=bioThfilt_norm;


    
Ncs=[ampfiltTh_norm,ampfiltAb_norm,ampfiltTh_no_norm,ampfiltAb_no_norm];
Bio=[bioThfilt_norm,bioAbfilt_norm];

NcsNum=4;BioNum=1;
 ampfilt1=Ncs(:,NcsNum);
 ampfilt2=Bio(:,BioNum);
NcsName={'Th wear','Ab wear','Th notch','Ab notch'};
BioName={'Th','Ab'};
NcsData=NcsName(NcsNum);
BioData=BioName(BioNum);
%%

opts3.tWinBR = 20; % Window on which br is estimated
opts3.tWin = 2; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.2;
    opts3.tWinVar=20;
    opts3.tWinpp=20;


 ampfilt2=ampfiltTh_norm;
tInd_normal=[1 fsDS*30+1;fsDS*60 fsDS*90;fsDS*120 fsDS*150]';
tInd_str=[fsDS*30 fsDS*60;fsDS*90 fsDS*120;fsDS*150 NaN]';
%tInd_str=[fsDS*90 fsDS*120;fsDS*90 fsDS*120;fsDS*150 NaN]';
ampfilt2_normal=[ampfilt2(tInd_normal(1,1):tInd_normal(2,1)),...
    ampfilt2(tInd_normal(1,2):tInd_normal(2,2)),ampfilt2(tInd_normal(1,3):tInd_normal(2,3))];
ampfilt2_str=[ampfilt2(tInd_str(1,1):tInd_str(2,1)),...
    ampfilt2(tInd_str(1,2):tInd_str(2,2)),ampfilt2(length(ampfilt2)-fsDS*30:end)];

[pk_nor,varFeature_nor] = brEst_Pk(ampfilt2_normal,fsDS,opts3);
[pk_str,varFeature_str] = brEst_Pk(ampfilt2_str,fsDS,opts3);

% var pp, BR, in-out time all normalized and *100   coefficent of Var (percent)
h(2)=figure;
subplot(3,1,1)
plot(tOff,brNcs,'LineWidth',1,'color','red');
hold on
plot(tOff,brBio,'LineWidth',1,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Breath Rate','FontSize',sz)
legend('NCS','chest belt','FontSize',sz)
title(fileName,'FontSize',sz)
ylim([0 80])

subplot(3,1,2)
plot(tOff,varBio(:,1),'color','red','LineWidth',1)
xlabel('time (s)','FontSize',sz)
ylabel('BR variation','FontSize',sz)
subplot(3,1,3)
plot(tOff,varBio(:,2),'color','green','LineWidth',1)
xlabel('time (s)','FontSize',sz)
ylabel('PP variation','FontSize',sz)

%%
%STFT 
% opt6.twin=fs*8;
% 
% fs=fsDS;
% [s,f,t] = stft(bioAbfilt_norm,fs,'window',hamming(fs*30,'periodic'),'FFTLength',fs*100);
% 
% s=abs(s(((length(f)/2+1):(length(f)-1)),:));
% f=f((length(f)/2+1):(length(f)-1));
% [M,I] =max(s);
% freq=f(I)*60;
% figure()
% plot(t,freq)
% 
% opt6.twinVar=4;
% Var_stft=zeros(length(t),1);
% for i=opt6.twinVar:length(t)
%     temp=zeros(opt6.twinVar,1);
%     for j= 1:opt6.twinVar
%     temp(j)=freq(i-j+1);
%     end
%     Var_stft(i)=var(temp);
% end
% Var_mean_stft=mean(Var_stft(Var_stft~=0));
% figure()
% plot(t,Var_stft)

%%
figName1 = [dataPath,'fig\',fileName,'waveform'];
print(h(1),[figName1,'.tiff'],'-dtiff','-r300');
savefig(h(1),[figName1,'.fig']);
figName2 = [dataPath,'fig\',fileName,'BRvar'];
print(h(2),[figName2,'.tiff'],'-dtiff','-r300');
savefig(h(2),[figName2,'.fig']);



BR=struct('brBio',brBio,'BioPk',BioPk,'varBio',varBio,'varBioFeature',varBioFeature,'brNcs',brNcs,'NcsPk',NcsPk,'varNcs',varNcs,'varNcsFeature',varNcsFeature,'NcsSelect',NcsData,'BioSelect',BioData,'period',period);
save([dataPath,'mfile\BRfeature\',fileName,'.mat'],'BR');
eval([fileName,'=','BR',';']);
