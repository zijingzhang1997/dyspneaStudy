clear all
dataPath=['C:\Sleep test\dyspnea\data\case1-6\'];
%x = categorical({'VarBRmean','Varppmean','VarInmean','VarExmean', 'meanBR', 'meanIn', 'meanEx' ,'BRcor',...
%   'ppcor','Incor','Excor'})';
delta=[];
delta_norm=[];
dysp=[];
fileInd={};
for i=['1','2','3','4','5','6']
    feat=[];
    fileTemp={};
    for j=['1','3']
        fileName=['Case',i,'Routine',j];
        [feat((end+1),:)]=loadData(dataPath,fileName);
        %[var((end+1),:),cor((end+1),:,:)]=loadData(dataPath,fileName);
        fileTemp{end+1}=fileName;
    end

    delta(end+1,:)=feat(1,:)-feat(2,:);
    delta_norm(end+1,:)=(feat(1,:)-feat(2,:))./feat(1,:).*100;
    dysp(end+1,:)=feat(2,:);
    fileInd{end+1,1}=fileTemp{1};
    fileInd{end,2}=fileTemp{2};
    
end

for i=['2','3','4','5','6']
    feat=[];
    fileTemp={};
    for j=['1','5']
        fileName=['Case',i,'Routine',j];
        [feat((end+1),:)]=loadData(dataPath,fileName);
        %[var((end+1),:),cor((end+1),:,:)]=loadData(dataPath,fileName);
        fileTemp{end+1}=fileName;
    end

    delta(end+1,:)=feat(1,:)-feat(2,:);
    delta_norm(end+1,:)=(feat(1,:)-feat(2,:))./feat(1,:).*100;
    dysp(end+1,:)=feat(2,:);
    fileInd{end+1,1}=fileTemp{1};
    fileInd{end,2}=fileTemp{2};
end

for i=['1','2','3','6']
     feat=[];
    fileTemp={};
    for j=['2','4']
        fileName=['Case',i,'Routine',j];
        [feat((end+1),:)]=loadData(dataPath,fileName);
        %[var((end+1),:),cor((end+1),:,:)]=loadData(dataPath,fileName);
        fileTemp{end+1}=fileName;
    end

    delta(end+1,:)=feat(1,:)-feat(2,:);
    delta_norm(end+1,:)=(feat(1,:)-feat(2,:))./feat(1,:).*100;
    dysp(end+1,:)=feat(2,:);
    fileInd{end+1,1}=fileTemp{1};
    fileInd{end,2}=fileTemp{2};
end


scale0=[7 3 4 4 2 5 5 4 1 4 5 7 4 3 4];
scale=[3 3 1 0 2 2 3 2 2 2 2 3 2 1 2];

%save(('C:\Sleep test\dyspnea\data\case1-6\data.mat'),'delta','delta_norm','scale','scale0');
% sum_delta=sum([delta(:,1:2) delta(:,8:9)],2);
% mean_delta_norm=mean([delta_norm(:,1:2) delta_norm(:,8:9)],2);

function [feature]=loadData(dataPath,fileName)
load([dataPath,fileName,'.mat']);
data=BR;


var=[data.varBioFeature.VarBR_mean 
    data.varBioFeature.Varpp_mean 
    data.varBioFeature.VarIn_mean 
    data.varBioFeature.VarEx_mean 
    data.varBioFeature.meanBR 
    data.varBioFeature.meanIn 
    data.varBioFeature.meanEx];

cor=data.correlationBio;
cor1=cor(:,1);
cor3=cor(:,3);
hrv=data.HRVfeature';
feature=[var; cor1;cor3; hrv];

end