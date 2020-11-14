clear all
dataPath=['C:\Sleep test\dyspnea\data\10_12\mfile\BRfeature\'];
%x = categorical({'VarBRmean','Varppmean','VarInmean','VarExmean', 'meanBR', 'meanIn', 'meanEx' ,'BRcor',...
%   'ppcor','Incor','Excor'})';
delta=[];
delta_norm=[];
fileInd={};
for i=[7,8,9,10,11]
    feat=[];
    fileTemp={};
    for j=['1','3']
        fileName=['Case',num2str(i),'Routine',j];
        [feat((end+1),:)]=loadDataNcs(dataPath,fileName);
        %[var((end+1),:),cor((end+1),:,:)]=loadData(dataPath,fileName);
        fileTemp{end+1}=fileName;
    end

    delta(end+1,:)=feat(1,:)-feat(2,:);
    delta_norm(end+1,:)=(feat(1,:)-feat(2,:))./feat(1,:).*100;
    fileInd{end+1,1}=fileTemp{1};
    fileInd{end,2}=fileTemp{2};
end
for i=[7,8,9,10,11]
    feat=[];
    fileTemp={};
    for j=['1','5']
        fileName=['Case',num2str(i),'Routine',j];
        [feat((end+1),:)]=loadDataNcs(dataPath,fileName);
        %[var((end+1),:),cor((end+1),:,:)]=loadData(dataPath,fileName);
        fileTemp{end+1}=fileName;
    end

    delta(end+1,:)=feat(1,:)-feat(2,:);
    delta_norm(end+1,:)=(feat(1,:)-feat(2,:))./feat(1,:).*100;
    fileInd{end+1,1}=fileTemp{1};
    fileInd{end,2}=fileTemp{2};
end
for i=[7,8,9,11]
    feat=[];
    fileTemp={};
    for j=['2','4']
        fileName=['Case',num2str(i),'Routine',j];
        [feat((end+1),:)]=loadDataNcs(dataPath,fileName);
        %[var((end+1),:),cor((end+1),:,:)]=loadData(dataPath,fileName);
        fileTemp{end+1}=fileName;
    end

    delta(end+1,:)=feat(1,:)-feat(2,:);
    delta_norm(end+1,:)=(feat(1,:)-feat(2,:))./feat(1,:).*100;
    fileInd{end+1,1}=fileTemp{1};
    fileInd{end,2}=fileTemp{2};
end
delta1=delta;
delta_norm1=delta_norm;
% for i=['2','3','4','5','6']
%     var=[];cor=[];
%     fileTemp={};
%     for j=['1','5']
%         fileName=['Case',i,'Routine',j];
%         [var((end+1),:),cor((end+1),:)]=loadData(dataPath,fileName);
%         %[var((end+1),:),cor((end+1),:,:)]=loadData(dataPath,fileName);
%         fileTemp{end+1}=fileName;
%     end
%     deltaVar=var(1,:)-var(2,:);
%     deltaVar_norm=(var(1,:)-var(2,:))./var(1,:).*100;
%     %deltaCor3=cor(1,:,3)-cor(2,:,3);
%     deltaCor3=cor(1,:)-cor(2,:);
%     deltaCor3_norm=(cor(1,:)-cor(2,:))./cor(1,:).*100;
%     deltaCor1=cor(1,:,1)-cor(2,:,1);
%     delta(end+1,:)=[deltaVar deltaCor3];
%     delta_norm(end+1,:)=[deltaVar_norm deltaCor3_norm];
%     fileInd{end+1,1}=fileTemp{1};
%     fileInd{end,2}=fileTemp{2};
% end
% 
% for i=['1','2','3','6']
%     var=[];cor=[];
%     fileTemp={};
%     for j=['2','4']
%         fileName=['Case',i,'Routine',j];
%         [var((end+1),:),cor((end+1),:)]=loadData(dataPath,fileName);
%         %[var((end+1),:),cor((end+1),:,:)]=loadData(dataPath,fileName);
%         fileTemp{end+1}=fileName;
%     end
%     deltaVar=var(1,:)-var(2,:);
%     deltaVar_norm=(var(1,:)-var(2,:))./var(1,:).*100;
%     %deltaCor3=cor(1,:,3)-cor(2,:,3);
%     deltaCor3=cor(1,:)-cor(2,:);
%     deltaCor3_norm=(cor(1,:)-cor(2,:))./cor(1,:).*100;
%     deltaCor1=cor(1,:,1)-cor(2,:,1);
%     delta(end+1,:)=[deltaVar deltaCor3];
%     delta_norm(end+1,:)=[deltaVar_norm deltaCor3_norm];
%     fileInd{end+1,1}=fileTemp{1};
%     fileInd{end,2}=fileTemp{2};
% end
% 
% sum_delta=sum([delta(:,1:4) delta(:,8:end)],2);
% mean_delta_norm=mean([delta_norm(:,1:4) delta_norm(:,8:end)],2);
% mean_delta_norm1=mean( delta_norm(:,8:11),2);
% mean_delta_norm2=mean(delta_norm(:,12:15),2);
% mean_delta_norm3=mean(delta_norm(:,16:19),2);
 scale0=[7 3 4 4 2 5 5 4 1 4 5 7 4 3 4]; %subject
 scale=[3 3 1 0 2 2 3 2 2 2 2 3 2 1 2];   %nurse 

scale2=[3,3,3,3,2,3,3,2,1,1,3,3,3,2];
scale=[scale scale2];

%save(('C:\Sleep test\dyspnea\data\case1-6\data.mat'),'delta','delta_norm','fileInd','sum_delta','mean_delta_norm','scale','scale0');
featureName=['VarBR','Varpp','VarIn','VarEx', 'meanBR',...
             'meanIn', 'meanEx' ,'BRcor1','ppcor1','Incor1','Excor1','BRcor3','ppcor3','Incor3','Excor3',...
            'meanHR','sdnn','rmsrr'];
% load(('C:\Sleep test\dyspnea\data\case1-6\data.mat'),'delta','delta_norm');
% delta=[delta ;delta1];
% delta_norm=[delta_norm; delta_norm1];
%save(('C:\Sleep test\dyspnea\data\case1-6\data.mat'),'delta','delta_norm','fileInd','scale','featureName');

function [feature]=loadDataBio(dataPath,fileName)
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

hrv=data.HRVfeature';
feature=[var;cor(:,1); cor(:,3); hrv];

end

function [feature]=loadDataNcs(dataPath,fileName)
load([dataPath,fileName,'.mat']);
data=BR;


var=[data.varNcsFeature.VarBR_mean 
    data.varNcsFeature.Varpp_mean 
    data.varNcsFeature.VarIn_mean 
    data.varNcsFeature.VarEx_mean 
    data.varNcsFeature.meanBR 
    data.varNcsFeature.meanIn 
    data.varNcsFeature.meanEx];

cor=data.correlationNcs;

hrv=data.HRVfeature';
feature=[var;cor(:,1); cor(:,3); hrv];

end