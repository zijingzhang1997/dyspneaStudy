%NCS 4 channels Bio 2 channels 
%select NCS by correlation to bio
%select bio by lower variation




function [NcsNum,BioNum,signNcs] = SelectData(Ncs,Bio,opts,fsDS)
[r1,p1] = corrcoef(Bio(:,2),Ncs(:,1)); 
r1 = r1(1,2);
[r2,p2] = corrcoef(Bio(:,2),Ncs(:,2)); 
r2 = r2(1,2); 
[r3,p3] = corrcoef(Bio(:,2),Ncs(:,3)); 
r3 = r3(1,2); 
[r4,p4] = corrcoef(Bio(:,2),Ncs(:,4));
r4 = r4(1,2); 
signNcs = sign([r1, r2, r3, r4]);
[rmax,idxCorr] = max(abs([r1,r2,r3,r4]));
NcsNum=idxCorr;

set(0,'DefaultFigureVisible','off')
% [~,~,~,varFeature1] = brEstAvg(Ncs(:,1),fsDS,opts);
% [~,~,~,varFeature2] = brEstAvg(Ncs(:,2),fsDS,opts);
% [~,~,~,varFeature3] = brEstAvg(Ncs(:,3),fsDS,opts);
% [~,~,~,varFeature4] = brEstAvg(Ncs(:,4),fsDS,opts);
% 
% 
% 
% [~,idxVar] = min(([sum(varFeature1(1:4)),sum(varFeature2(1:4)),sum(varFeature3(1:4)),sum(varFeature4(1:4))]));
%NcsNum=idxVar;

NcsName={'Th wear','Ab wear','Th notch','Ab notch'};
BioName={'Th','Ab'};

fprintf(' NCS : use %s ',NcsName{NcsNum});



[brbioTh,bioThPk,varbioTh,varbioThFeature] = brEstAvg(Bio(:,1),fsDS,opts);
[brbioAb,bioAbPk,varbioAb,varbioAbFeature] = brEstAvg(Bio(:,2),fsDS,opts);

[~,idx] = min(([sum(varbioThFeature(1:4)),sum(varbioAbFeature(1:4))]));
BioNum=idx;
fprintf(' Bio : use %s \n',BioName{BioNum});
set(0,'DefaultFigureVisible','on')

end