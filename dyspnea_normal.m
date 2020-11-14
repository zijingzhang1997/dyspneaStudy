
dataPath=['C:\Sleep test\dyspnea\data\all\'];

% for i=['4','5']
%     for j=['1','2','3','4','5','6']
%         fileName=['Case',i,'Routine',j]
%         dyspnea_cor(dataPath,fileName);
%     end
% end
for i=[2:15]
    for j=[1]
        CaseNum=[i,j];
        dyspnea_nor(dataPath,CaseNum,1);
        dyspnea_nor(dataPath,CaseNum,2);
    end
end






function dyspnea_nor(dataPath,CaseNum,opt)
close all
% dataPath=['C:\Sleep test\dyspnea\data\10_6\'];
% fileName='Case3Routine4';
fileName=['Case',num2str(CaseNum(1)),'Routine',num2str(CaseNum(2))];
fs=10e3;
fsDS=500;

%toff=[5:295]';

tNoise=[];


fsbio=1e3;
fsDSbio=500;
filePathName = [dataPath,fileName,'.tdms'];
filePathName_m = [dataPath,fileName,'.mat'];
if ~exist(filePathName_m,'file') & ~exist(filePathName,'file')
    return
end
if ~exist(filePathName_m,'file')
   convertTDMS(true,filePathName);
end
load(filePathName_m);

% if exist([dataPath,'mfile_old\BRfeature\',fileName,'.mat']) ~=0
%     load([dataPath,'mfile_old\BRfeature\',fileName,'.mat'],'BR');
%     fprintf('substitue time');
%     toff=[BR.period(1):BR.period(2)]';
% end

load([dataPath,fileName,'.mat']);
ampTh=ConvertedData.Data.MeasuredData(3).Data;%  rx1 thorax 
bioTh=ConvertedData.Data.MeasuredData(20).Data;  %ch2 biopac thorax 
ampTh_no=ConvertedData.Data.MeasuredData(13).Data;  %13  tx3
ampAb_no=ConvertedData.Data.MeasuredData(18).Data;   %18 tx4
ampAb=ConvertedData.Data.MeasuredData(8).Data; %rx2 abdomen
bioAb=ConvertedData.Data.MeasuredData(21).Data;  %ch3 biopac abdomen
PPG=ConvertedData.Data.MeasuredData(22).Data;  %19 biopac ch1
ECG=ConvertedData.Data.MeasuredData(19).Data;

len=round(length(bioAb)/fsbio);
if opt==1
toff=[10:round((len-15)/2)+10]';
end
 if opt==2
toff=[round((len-15)/2)+10:len-5]';
end

period=[toff(1),toff(end)];
[HRV] = HRVecg(dataPath,fileName,toff);
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





%%  breath rate estimation  peak-to-peak

%ampfilt1=ampfiltTh_norm;
%ampfilt2=bioThfilt_norm;

opts3.tWinBR = 20; % Window on which br is estimated
opts3.tWin = 2; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.2;
    opts3.tWinVar=20;
    opts3.tWinpp=20;
    
Ncs=[ampfiltTh_norm,ampfiltAb_norm,ampfiltTh_no_norm,ampfiltAb_no_norm];
Bio=[bioThfilt_norm,bioAbfilt_norm];
[NcsNum,BioNum,signNcs] = SelectData(Ncs,Bio,opts3,fsDS);

 ampfilt1=Ncs(:,NcsNum)*signNcs(NcsNum);
 ampfilt2=Bio(:,BioNum);

NcsName={'Th wear','Ab wear','Th notch','Ab notch'};
BioName={'Th','Ab'};
NcsData=NcsName(NcsNum);
BioData=BioName(BioNum);

h(1)=figure;
sz=10;
titleName=[fileName,'  Ncs:',char(NcsData),'  Bio:',char(BioData)];

subplot(6,1,1);
plot(tOff,ampfiltTh_norm,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
%xlim([0 220])


title('NCS Th','FontSize',sz);

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

sgtitle(titleName)
set(gcf,'Position',[100,200,1000,800]);
figName = [dataPath,'fig\',fileName,'waveform'];
%%

opts3.tWinBR = 20; % Window on which br is estimated
opts3.tWin = 2; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.2;
    opts3.tWinVar=20;
    opts3.tWinpp=20;

[cor1,h(2),ampfilt1New,temp1,SD1] = brEstCor(ampfilt1,fsDS,CaseNum);
[cor2,h(3),ampfilt2New,temp2,SD2] = brEstCor(ampfilt2,fsDS,CaseNum);

[brNcs,NcsPk,varNcs,varNcsFeature,h(4)] = brEstAvg(ampfilt1New,fsDS,opts3);


[brBio,BioPk,varBio,varBioFeature,h(5)] = brEstAvg(ampfilt2New,fsDS,opts3);

tOffNewNcs=((0:(length(ampfilt1New(:,1))-1))/fsDSbio)';
tOffNewBio=((0:(length(ampfilt2New(:,1))-1))/fsDSbio)';
[rNCSBio,] = corrcoef(ampfilt1,ampfilt2);
rNCSBio = rNCSBio(1,2); 

% var pp, BR, in-out time all normalized and *100   coefficent of Var (percent)

% subplot(3,1,1)
% plot(tOffNewNcs,brNcs,'LineWidth',1,'color','red');
% hold on
% plot(tOffNewBio,brBio,'LineWidth',1,'color','green');
% xlabel('time (s)','FontSize',sz)
% ylabel('Breath Rate','FontSize',sz)
% legend('NCS','chest belt','FontSize',sz)
% title(fileName,'FontSize',sz)
% ylim([0 80])
% 
% subplot(3,1,2)
% plot(tOffNewBio,varBio(:,1),'color','red','LineWidth',1)
% xlabel('time (s)','FontSize',sz)
% ylabel('BR variation','FontSize',sz)
% subplot(3,1,3)
% plot(tOffNewBio,varBio(:,2),'color','green','LineWidth',1)
% xlabel('time (s)','FontSize',sz)
% ylabel('PP variation','FontSize',sz)
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
fileName=[fileName,'_',num2str(opt)];
figName1 = [dataPath,'fig_normal\',fileName,'waveform'];
print(h(1),[figName1,'.tiff'],'-dtiff','-r300');
savefig(h(1),[figName1,'.fig']);
figName2 = [dataPath,'fig_normal\',fileName,'corrNcs'];
print(h(2),[figName2,'.tiff'],'-dtiff','-r300');
savefig(h(2),[figName2,'.fig']);
figName3 = [dataPath,'fig_normal\',fileName,'corrBio'];
print(h(3),[figName3,'.tiff'],'-dtiff','-r300');
savefig(h(3),[figName3,'.fig']);
figName4 = [dataPath,'fig_normal\',fileName,'peakNcs'];
print(h(4),[figName4,'.tiff'],'-dtiff','-r300');
savefig(h(4),[figName4,'.fig']);
figName5 = [dataPath,'fig_normal\',fileName,'peakBio'];
print(h(5),[figName5,'.tiff'],'-dtiff','-r300');
savefig(h(5),[figName5,'.fig']);

BR=struct('brBio',brBio,'BioPk',BioPk,'varBio',varBio,'varBioFeature',varBioFeature,'brNcs',brNcs,'NcsPk',NcsPk,...
    'varNcs',varNcs,'varNcsFeature',varNcsFeature,'correlationNcs',cor1,'correlationBio',cor2,'NcsSelect',NcsData,...
    'BioSelect',BioData,'HRVfeature',HRV,'period',period,'outlier',[temp1;temp2],'rNCSBio',rNCSBio,'SD',[SD1 SD2]);
%varFeature : VarBR_mean VarPP_mean VarIn_mean VarEx_mean meanBR meanPP meanIn meanEx
%SD standard deviation(1:4) mean(5:8) for the whole routine
%rNCSBio  correlation of Bio and Ncs 
save([dataPath,'mfile_normal\BRfeature\',fileName,'.mat'],'BR');
eval([fileName,'=','BR',';']);
end