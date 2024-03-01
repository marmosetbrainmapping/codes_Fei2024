load('bilateralFCmap.mat');
load('/zfs/scratch/qihang/R2L_SC_strength.mat');
load('/zfs/scratch/qihang/cortex_mask.mat');
bilateralFCmap(logical(bilateralFCmap==1)) = nan;
bilateralFCmap(logical(cortex_mask==0)) = nan;


%%
Allennames = ["FRP","MOp","MOs","SSp-n","SSp-bfd","SSp-ll","SSp-m","SSp-ul","SSp-tr","SSp-un","SSs","GU","VISC","AUDd","AUDp","AUDpo","AUDv","VISal","VISam","VISl","VISp","VISpl","VISpm","VISli","VISpor","ACAd","ACAv","PL","ILA","ORBl","ORBm","ORBvl","AId","AIp","AIv","RSPagl","RSPd","RSPv","VISa","VISrl","TEa","PERI","ECT"];
%%
load('ccfv3_mask.mat');
%%
ccfv3_names = [];
for i = 1:length(ccfv3_mask)
    ccfv3_names = [ccfv3_names,ccfv3_mask(i).name];
end
%%
[commonnames,iallen,iccfv3] = intersect(Allennames,ccfv3_names);
%%
bilateralFC = zeros(length(commonnames),1);
bilateralSC = zeros(length(commonnames),1);

for i = 1:length(commonnames)
    
    ccfv3_idx = iccfv3(i);
    roi = double(ccfv3_mask(ccfv3_idx).mask);
    roi(logical(roi==0)) = nan;
    FC = roi.*bilateralFCmap;
    FC = mean(FC,[1,2],'omitnan');
    
    allen_idx = iallen(i);
    SC = connectivity(allen_idx,allen_idx);
    
    bilateralFC(i) = FC;
    bilateralSC(i) = SC;
    
end

%%
figure
hold on

x = bilateralSC;
y = bilateralFC;

f = fit(x,y,'poly1');
k = f.p1;
b = f.p2;
x_ = linspace(0,1,50);
y_ = k*x_ + b;
plot(x_,y_,'black','LineWidth',3);
s = scatter(x,y,80,'o','filled','MarkerEdgeColor',[1,0.4,0.4],'MarkerFaceColor',[1,0.4,0.4]);

[R,P,RL,RU] = corrcoef(x,y,'Alpha',0.05);
r = R(1,2);
pvalue = P(1,2);
title(['r=',num2str(round(r*1000)/1000),',p=',num2str(round(pvalue*1000)/1000)]);
ylabel('homotopic FC');
xlabel('homotopic SC');
set(gca,'Fontsize',15,'Color','none');
set(gcf,'Position',[0,0,550,500]);

xlim([0,1]);
ylim([0,1]);



