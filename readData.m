clear all
dataPath=['C:\Sleep test\dyspnea\data\all\mfile\BRfeature\'];
%x = categorical({'VarBRmean','Varppmean','VarInmean','VarExmean', 'meanBR', 'meanIn', 'meanEx' ,'BRcor',...
%   'ppcor','Incor','Excor'})';
deltaBio=[];
deltaBio_norm=[];
fileIndBio={};
deltaNcs=[];
deltaNcs_norm=[];
fileIndNcs={};
diff=[];
for i=[1:15]
    featNcs=[];
    featBio=[];
    fileTemp={};
    for j=[1,3]
        fileName=['Case',num2str(i),'Routine',num2str(j)];
        [featNcs((end+1),:)]=loadDataNcs(dataPath,fileName);
        [featBio((end+1),:)]=loadDataBio(dataPath,fileName);
        diff(end+1)=sum(featNcs((end),1:4)-featBio((end),1:4));
        fileTemp{end+1}=fileName;
    end
    
    deltaBio(end+1,:)=featNcs(1,:)-featNcs(2,:);
    deltaBio_norm(end+1,:)=(featNcs(1,:)-featNcs(2,:))./featNcs(1,:).*100;
    fileIndBio{end+1,1}=fileTemp{1};
    fileIndBio{end,2}=fileTemp{2};
    
 

    
    deltaNcs(end+1,:)=featBio(1,:)-featBio(2,:);
    deltaNcs_norm(end+1,:)=(featBio(1,:)-featBio(2,:))./featBio(1,:).*100;

    


  
end

for i=[2:15]
    featNcs=[];
    featBio=[];
    fileTemp={};
    for j=[1,5]
        fileName=['Case',num2str(i),'Routine',num2str(j)];
        [featNcs((end+1),:)]=loadDataNcs(dataPath,fileName);
        [featBio((end+1),:)]=loadDataBio(dataPath,fileName);
        diff(end+1)=sum(featNcs((end),1:4)-featBio((end),1:4));
        fileTemp{end+1}=fileName;
    end
    
    deltaBio(end+1,:)=featNcs(1,:)-featNcs(2,:);
    deltaBio_norm(end+1,:)=(featNcs(1,:)-featNcs(2,:))./featNcs(1,:).*100;
    fileIndBio{end+1,1}=fileTemp{1};
    fileIndBio{end,2}=fileTemp{2};

    
    deltaNcs(end+1,:)=featBio(1,:)-featBio(2,:);
    deltaNcs_norm(end+1,:)=(featBio(1,:)-featBio(2,:))./featBio(1,:).*100;


  
end

%normal breathing half 1 2
dataPathNorm=['C:\Sleep test\dyspnea\data\all\mfile_normal\BRfeature\'];

% for i=[1:15]
%     featNcs=[];
%     featBio=[];
%     fileTemp={};
%     for j=[1,2]
%         fileName=['Case',num2str(i),'Routine1_',num2str(j)];
%         [featNcs((end+1),:)]=loadDataNcs(dataPathNorm,fileName);
%         [featBio((end+1),:)]=loadDataBio(dataPathNorm,fileName);
%         diff(end+1)=sum(featNcs((end),1:4)-featBio((end),1:4));
%         fileTemp{end+1}=fileName;
%     end
%     
%     deltaBio(end+1,:)=featNcs(1,:)-featNcs(2,:);
%     deltaBio_norm(end+1,:)=(featNcs(1,:)-featNcs(2,:))./featNcs(1,:).*100;
%     fileIndBio{end+1,1}=fileTemp{1};
%     fileIndBio{end,2}=fileTemp{2};
%     
%  
% 
%     
%     deltaNcs(end+1,:)=featBio(1,:)-featBio(2,:);
%     deltaNcs_norm(end+1,:)=(featBio(1,:)-featBio(2,:))./featBio(1,:).*100;
% 
%    
% end

 
 scale13=[3 1 2 2 1 2 3 3 3 3 2 2 2 3 2];
 %scale13_n=[2 3 1 1 2 2 3 3 1 1 1 1 2 3 2];
 scale13_n=[2 3 1 1 2 1 2 3 1 1 1 2 2 3 1];
%subject
 scale15=[2 2 1 2 2 3 3 2 1 1 1 1 1 1 ];   %nurse 
 scale15_n=[3 2 1 2 2 1 3 1 1 1 2 1 1 1];
 scale24=[];
scale11=zeros(1,15);
scale_n=[scale13_n scale15_n];
scale=[scale13 scale15];
% scale_n=[scale13_n scale15_n scale11];
% scale=[scale13 scale15 scale11];
% std_br=[];
% std_pp=[];
% std_in=[];
% std_ex=[];
% for i=[1:15]
%     featBio=[];
%     fileTemp={};
%     br=[];
%     pp=[];
%     in=[];
%     ex=[];
%     for j=[1:5]
%         fileName=['Case',num2str(i),'Routine',num2str(j)];
%    
%         [featBio((end+1),:)]=loadDataBio(dataPath,fileName);
%         br(j,1)=featBio(end,19);br(j,2)=featBio(end,5);
%         pp(j,1)=featBio(end,20);pp(j,2)=featBio(end,20);
%         in(j,1)=featBio(end,21);in(j,2)=featBio(end,6);
%         ex(j,1)=featBio(end,22);ex(j,2)=featBio(end,7);
%         fileTemp{end+1}=fileName;
%     end
%     
%     std_br(end+1,:,:)=br(:,:);
%     std_pp(end+1,:,:)=pp(:,:);
%     std_in(end+1,:,:)=in(:,:);
%     std_ex(end+1,:,:)=ex(:,:);
% end
% deltaBio=deltaBio(1:29,:);
% deltaBio_norm=deltaBio_norm(1:29,:);
% deltaNcs=deltaNcs(1:29,:);
% deltaNcs_norm=deltaNcs_norm(1:29,:);
% fileIndBio=fileIndBio(1:29,:);


%featureName=['VarBR','Varpp','VarIn','VarEx', 'meanBR',...
%              'meanIn', 'meanEx' ,'BRcor1','ppcor1','Incor1','Excor1','BRcor3','ppcor3','Incor3','Excor3',...
%          'meanHR','sdnn','rmsrr','stdBR','stdpp','stdIn','stdEx'];

sum_delta(:,1)=sum(deltaNcs(:,1:4),2);
sum_delta(:,2)=sum(deltaNcs(:,12:15),2);
scale=scale';
 delta2=[deltaNcs ;deltaBio];
 scale2=[scale scale];
 scale2_n=[scale_n scale_n];
 
%save(('C:\Sleep test\dyspnea\data\all\data_norm.mat'),'deltaBio','deltaNcs','delta2','fileIndBio','scale2_n','scale2');
save(('C:\Sleep test\dyspnea\data\all\data_wo.mat'),'deltaBio','deltaNcs','delta2','fileIndBio','scale2_n','scale2');
 
 
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
std=data.SD(1:4,2);
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
std=data.SD(1:4,1);
hrv=data.HRVfeature';
feature=[var;cor(:,1); cor(:,3); hrv];

end