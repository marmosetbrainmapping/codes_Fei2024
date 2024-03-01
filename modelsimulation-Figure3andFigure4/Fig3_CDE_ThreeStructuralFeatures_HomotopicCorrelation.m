
%%
lines = ["-","--",":","-.","-","--"];
markers = ["o","^","diamond","*",">","<"];
% colors = ["r","g","b","c","y","m"];
colors = {[1,0,0],[114,113,113]/255,[0,0,1]};
%%
lines = ["-","--",":","-.","-","--"];
markers = ["o","^","diamond","*",">","<"];

heterocolor = [231,143,57]/255;
homocolor = [46,167,224]/255;
pcontracolor = [45,92,156]/255;
sz = 10;
%%
figure
hold on
load('HomoParams_Nregion2.mat','Cmeans','noises','overlaphomo_range');
legends = [];
for i = 1:length(noises)-1
    x = overlaphomo_range;
    y = mean(abs(Cmeans(:,i,:)),1);
    plot(x,squeeze(y),'Color',colors{i},'LineWidth',1.5,'LineStyle',lines(i),'Marker',markers(i),'MarkerSize',sz,'MarkerEdgeColor',colors{i},'MarkerFaceColor','none');
    legends = [legends,string(num2str(noises(i)))];
end
ylim([0.45,0.63]);
yticks(0.4:0.1:0.8);
xticks(0:0.2:1);
set(gca,'Fontsize',15);
set(gca,'Color','none');
l = legend(legends);
set(l,'Position',[0.7185    0.1437    0.1643    0.1690]);
set(l,'Color','none');
set(gcf,'Position',[0,0,500,400]);

%%
figure
hold on
load('HeteroParams_Nregion2.mat','Cmeans','noises','overlaphetero_range');
legends = [];
for i = 1:length(noises)-1
    x = 1-overlaphetero_range;
    y = mean(abs(Cmeans(:,i,:)),1);
    plot(x,squeeze(y),'Color',colors{i},'LineWidth',1.5,'LineStyle',lines(i),'Marker',markers(i),'MarkerSize',sz,'MarkerEdgeColor',colors{i},'MarkerFaceColor','none');
    legends = [legends,string(num2str(noises(i)))];
end
ylim([0.45,0.63]);
yticks(0.4:0.1:0.8);
xticks(0:0.2:1);
set(gca,'Fontsize',15);
set(gca,'Color','none');
l = legend(legends);
% set(l,'Position',[0.7185    0.1437    0.1643    0.1690]);
set(l,'Color','none');
set(gcf,'Position',[0,0,500,400]);

%%
figure
hold on
load('PcontrasParams_Nregion2.mat','Cmeans','noises','pcontras_range');
legends = [];
for i = 1:length(noises)-1
    x = pcontras_range;
    y = mean(abs(Cmeans(:,i,:)),1);
    plot(x,squeeze(y),'Color',colors{i},'LineWidth',1.5,'LineStyle',lines(i),'Marker',markers(i),'MarkerSize',sz,'MarkerEdgeColor',colors{i},'MarkerFaceColor','none');    
    legends = [legends,string(num2str(noises(i)))];
end
ylim([0.45,0.63]);
yticks(0.5:0.1:0.8);
xticks(0:0.2:1);
set(gca,'Fontsize',15);
set(gca,'Color','none');
l = legend(legends);
set(l,'Position',[0.7185    0.1437    0.1643    0.1690]);
set(l,'Color','none');
set(gcf,'Position',[0,0,500,400]);

%%

