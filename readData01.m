clear all
dataPath=['C:\Sleep test\dyspnea\data\all\mfile\BRfeature\'];
%x = categorical({'VarBRmean','Varppmean','VarInmean','VarExmean', 'meanBR', 'meanIn', 'meanEx' ,'BRcor',...
%   'ppcor','Incor','Excor'})';
 featNcs=[];
 featBio=[];
 fileInd=[];
for i=[1:15]
  
    for j=[1,3,5]
        fileName=['Case',num2str(i),'Routine',num2str(j)];
        [featNcs((end+1),:)]=loadDataNcs(dataPath,fileName);
        [featBio((end+1),:)]=loadDataBio(dataPath,fileName);
        
        fileInd{end+1}=fileName;
    end
    
    
end

featBio(3,:)=[];
featNcs(3,:)=[];
fileInd(3)=[];
 
 scale3_n=[2 3 1 1 2 1 2 3 1 1 1 2 2 3 1];
 scale3=[3 1 2 2 1 2 3 3 3 3 2 2 2 3 2];
 scale5=[0 2 2 1 2 2 3 3 2 1 1 1 1 1 1];
 scale5_n=[0 3 2 1 2 2 1 3 1 1 1 2 1 1 1];
 scale0=zeros(1,15);
 
 scale=[scale0; scale3;scale5];
 scale_n=[scale0; scale3_n;scale5_n];
 scale=reshape(scale,[1,45]);
 scale_n=reshape(scale_n,[1,45]);
 scale_n(3)=[];
 scale(3)=[];
 feat2=[featNcs ;featBio];
 scale2=[scale scale];
 scale2_n=[scale_n scale_n];
save(('C:\Sleep test\dyspnea\data\all\data01.mat'),'featBio','featNcs','scale','feat2','scale2','scale2_n');
%featureName=['VarBR','Varpp','VarIn','VarEx', 'meanBR',...
%              'meanIn', 'meanEx' ,'BRcor1','ppcor1','Incor1','Excor1','BRcor3','ppcor3','Incor3','Excor3',...
%          'meanHR','sdnn','rmsrr','stdBR','stdpp','stdIn','stdEx'];


function [feature]=loadDataBio(dataPath,fileName)
if ~exist ([dataPath,fileName,'.mat'],'file')
    feature=zeros(18,1);
    return
end
load([dataPath,fileName,'.mat']);
data=BR;


var=data.varBioFeature;

cor=data.correlationBio;

hrv=data.HRVfeature';

feature=[var;cor(:,1); cor(:,3); hrv ];

end

function [feature]=loadDataNcs(dataPath,fileName)
if ~exist ([dataPath,fileName,'.mat'],'file')
    feature=zeros(18,1);
    return
end
load([dataPath,fileName,'.mat']);
data=BR;


var=data.varNcsFeature;

cor=data.correlationNcs;

hrv=data.HRVfeature';
feature=[var;cor(:,1); cor(:,3); hrv];

end