corr_thres = 0.85;
%%
load(['ratioParams_thres',num2str(corr_thres),'.mat']);
colormap(parula(length(noises)-1));
cmap = colormap;
cmap = flipud(cmap);

% cmap = colormap(cbrewer2('RdBu', length(noises)-1));
% cmap = flipud(cmap);

lines = ["-","--",":","-.","--",":"];
markers = ["o","*","diamond","^",">","<"];
sz = 8;
%% Fig. 4F Strength of correlation
figure
hold on
for i = 1:length(noises)-1

    x = ratios;
    y = mean(abs(Cmeans(:,i,:)),1);
    plot(x,squeeze(y)/max(y),'Color',cmap(i,:),'LineWidth',1.5,'MarkerEdgeColor',cmap(i,:),'Marker',markers(i),'MarkerSize',sz,'LineStyle',lines(i),'MarkerFaceColor','none');

end
set(gca,'Color','none');
set(gca,'FontSize',15);
set(gcf,'Position',[0,0,600,400]);

%% Fig. 4G Variability of correlation
figure
hold on
for i = 1:length(noises)-1

    x = ratios;
    y = mean(abs(Cstds(:,i,:)),1);
    plot(x,squeeze(y)/max(y),'Color',cmap(i,:),'LineWidth',1.5,'MarkerEdgeColor',cmap(i,:),'Marker',markers(i),'MarkerSize',sz,'LineStyle',lines(i),'MarkerFaceColor','none');

end
set(gca,'Color','none');
set(gca,'FontSize',15);
set(gcf,'Position',[0,0,600,400]);

%% Fig. 4E Duration of synchronous periods
figure
hold on
legends = [];
for i = 1:length(noises)-1

    x = ratios;
    y = mean(abs(psyncs(:,i,:)),1);
    plot(x,squeeze(y)/max(y),'Color',cmap(i,:),'LineWidth',1.5,'MarkerEdgeColor',cmap(i,:),'Marker',markers(i),'MarkerSize',sz,'LineStyle',lines(i),'MarkerFaceColor','none');
    legends = [legends,string(num2str(noises(i)))];

end
l = legend(legends);
set(l,'Color','none');
set(gca,'Color','none');
set(gca,'FontSize',15);
set(gcf,'Position',[0,0,600,400]);

%%


