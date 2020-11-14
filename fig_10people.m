S1 = load('C:\Sleep test\Data\624\final\edwin');
S2 = load('C:\Sleep test\Data\624\final\thomas');
S3 = load('C:\Sleep test\Data\624\final\tejas');
S4 = load('C:\Sleep test\Data\624\final\xiaonan');
S5 = load('C:\Sleep test\Data\624\final\yunlei');
S6 = load('C:\Sleep test\Data\624\final\pragya_data');
S7 = load('C:\Sleep test\Data\624\final\jianlin');
S8 = load('C:\Sleep test\Data\624\final\volunteer');
S9 = load('C:\Sleep test\Data\624\final\huangqingbo');
S10 = load('C:\Sleep test\Data\624\final\zijing');
% brBio=S.brBio;
% brBio_1=S.brBio_1;
% brNcs=S.brNcs;
% brNcs_1=S.brNcs_1;
fsDS=500;

brBio_1=[S1.brBio_1;S2.brBio_1;S3.brBio_1;S4.brBio_1;S5.brBio_1;S6.brBio_1;S7.brBio_1;S8.brBio_1;S9.brBio_1;S10.brBio_1];
brNcs_1=[S1.brNcs_1;S2.brNcs_1;S3.brNcs_1;S4.brNcs_1;S5.brNcs_1;S6.brNcs_1;S7.brNcs_1;S8.brNcs_1;S9.brNcs_1;S10.brNcs_1];
brBio=[S1.brBio,S2.brBio,S3.brBio,S4.brBio,S5.brBio,S6.brBio,S7.brBio,S8.brBio,S9.brBio,S10.brBio];
brNcs=[S1.brNcs,S2.brNcs,S3.brNcs,S4.brNcs,S5.brNcs,S6.brNcs,S7.brNcs,S8.brNcs,S9.brNcs,S10.brNcs];

%% 
pkLS = regress(brBio_1,[ones(length(brNcs_1),1) brNcs_1]);
pkBaMean = mean(brBio_1 - brNcs_1)
pkBaStd = std(brBio_1 - brNcs_1)
pkBaStdLim = [pkBaStd*1.96+pkBaMean, -pkBaStd*1.96+pkBaMean]; 
[rPkVol,pPkVol] = corrcoef(brNcs_1,brBio_1);
mean=(brNcs+brBio)./2;
diff=brNcs-brBio;

%% 


sz=10;a=[5.8,62.7,75,164,194.5,length(brNcs(:,1))/fsDS];
figure
subplot(1,2,1);
brNcs1=reshape(brNcs(a(1)*fsDS:a(2)*fsDS,:),[],1);brBio1=reshape(brBio(a(1)*fsDS:a(2)*fsDS,:),[],1);
brNcs2=reshape(brNcs(a(2)*fsDS:a(3)*fsDS,:),[],1);brBio2=reshape(brBio(a(2)*fsDS:a(3)*fsDS,:),[],1);
brNcs3=reshape(brNcs(a(3)*fsDS:a(4)*fsDS,:),[],1);brBio3=reshape(brBio(a(3)*fsDS:a(4)*fsDS,:),[],1);
brNcs4=reshape(brNcs(a(4)*fsDS:a(5)*fsDS,:),[],1);brBio4=reshape(brBio(a(4)*fsDS:a(5)*fsDS,:),[],1);
brNcs5=reshape(brNcs(a(5)*fsDS:a(6)*fsDS,:),[],1);brBio5=reshape(brBio(a(5)*fsDS:a(6)*fsDS,:),[],1);
scatter(brNcs1,brBio1,'o','filled'); % make into different parts
hold on
scatter(brNcs2,brBio2,'+');
scatter(brNcs3,brBio3,'*');
scatter(brNcs4,brBio4,'x');
scatter(brNcs5,brBio5,'o','filled');

xVol = 0:10:60;
plot(xVol,pkLS(1)+pkLS(2)*xVol,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
text(50,30,['\bf', 'r=',num2str(rPkVol(1,2),3)],'FontSize',11)
xlabel('NCS BR (BPM)','FontSize',sz);
ylabel('BIOPAC BR (BPM)','FontSize',sz);
ylim([0 60]);

subplot(1,2,2);
mean1=reshape(mean(a(1)*fsDS:a(2)*fsDS,:),[],1);diff1=reshape(diff(a(1)*fsDS:a(2)*fsDS,:),[],1);
mean2=reshape(mean(a(2)*fsDS:a(3)*fsDS,:),[],1);diff2=reshape(diff(a(2)*fsDS:a(3)*fsDS,:),[],1);
mean3=reshape(mean(a(3)*fsDS:a(4)*fsDS,:),[],1);diff3=reshape(diff(a(3)*fsDS:a(4)*fsDS,:),[],1);
mean4=reshape(mean(a(4)*fsDS:a(5)*fsDS,:),[],1);diff4=reshape(diff(a(4)*fsDS:a(5)*fsDS,:),[],1);
mean5=reshape(mean(a(5)*fsDS:a(6)*fsDS,:),[],1);diff5=reshape(diff(a(5)*fsDS:a(6)*fsDS,:),[],1);


scatter(mean1,diff1,'o','filled'); % Normal
xVol = 0:10:60;
hold on
scatter(mean2,diff2,'+');
scatter(mean3,diff3,'*');
scatter(mean4,diff4,'x');
scatter(mean5,diff5,'o','filled');
plot(xVol,pkBaMean.*ones(length(xVol),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
plot(xVol,pkBaStdLim(1).*ones(length(xVol),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
plot(xVol,pkBaStdLim(2).*ones(length(xVol),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);
xlabel('mean BR (BPM)','FontSize',sz);
ylabel('difference BR (BPM)','FontSize',sz);
ylim([-40 40]);xlim([0 60]);
legend('normal breath','hold breath','deep breath','fast breath','normal breath 2','NumColumns',3);

text(50,pkBaMean+5,['\bf', num2str(pkBaMean,3)],'FontSize',10)
text(50,pkBaStdLim(1)+5,['\bf', num2str(pkBaStdLim(1),3)],'FontSize',10)
text(50,pkBaStdLim(2)+5,['\bf', num2str(pkBaStdLim(2),3)],'FontSize',10)
