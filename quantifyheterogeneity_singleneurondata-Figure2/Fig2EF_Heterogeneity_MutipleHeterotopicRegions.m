load('ipsi_neuron_to_region_count.mat');
ipsi_count = counts;
load('contra_neuron_to_region_count.mat');
contra_count = counts;

load('ITcortexneuronname.mat');
load('cortexregionnames.mat');

PFCregions = unique(ITcortexneuronname);
Cortexregions = (cortexregionnames);
%
PFCinds = zeros(length(PFCregions),1);
for i = 1:length(PFCregions)
    PFCregion = PFCregions(i);
    ind = find(Cortexregions==PFCregion);
    PFCinds(i) = ind;
end
%
ipsi_count_pfc = ipsi_count(:,PFCinds);
contra_count_pfc = contra_count(:,PFCinds);
%
neuroninds = cell(length(PFCinds),1);
for i = 1:length(PFCinds)
    PFCregion = PFCregions(i);
    inds = find(ITcortexneuronname==PFCregion);
    neuroninds{i,1} = inds;
end
%
ipsi_count_neurons = cell(length(PFCinds),1);
contra_count_neurons = cell(length(PFCinds),1);
for i = 1:length(PFCinds)
    inds = neuroninds{i};
    ipsi_count_neurons{i,1} = ipsi_count(inds,:);
    contra_count_neurons{i,1} = contra_count(inds,:);
end


ipsi_proj_neurons = cell(length(PFCregions),1);
contra_proj_neurons = cell(length(PFCregions),1);
for i = 1:length(PFCregions)
    i_projs = ipsi_count_neurons{i};
    ipsi_proj_neurons{i} = double(i_projs>0);
    i_projs = contra_count_neurons{i};
    contra_proj_neurons{i} = double(i_projs>0);
end

ipsi_intersect_count = zeros(length(PFCregions),length(Cortexregions)-1);
contra_intersect_count = zeros(length(PFCregions),length(Cortexregions)-1);
for j = 1:size(ipsi_intersect_count,1)
    regionid = j;
    ipsi_data = ipsi_proj_neurons{regionid};
    ipsi_data(:,PFCinds(regionid)) = [];
    contra_data = contra_proj_neurons{regionid};
    contra_data(:,PFCinds(regionid)) = [];
    n_neuron = size(ipsi_data,1);
    n_ipsi = sum(sum(ipsi_data,2)>0);
    n_contra = sum(sum(contra_data,2)>0);
    for i = 1:size(ipsi_intersect_count,2)
        count = sum(sum(ipsi_data,2)==i);
        ipsi_intersect_count(j,i) = count/n_ipsi;
        count = sum(sum(contra_data,2)==i);
        contra_intersect_count(j,i) = count/n_contra;
    end
end

x1 = ipsi_intersect_count;
x2 = contra_intersect_count;
for i = 1:size(x1,1)
    x1(i,:) = x1(i,:)/sum(x1(i,:));
    x2(i,:) = x2(i,:)/sum(x2(i,:));
end


%% Fig. 2E
figure
hold on

plot(nanmedian(x1(:,:),1)/nansum(nanmedian(x1(:,:),1))*100,'Color','r','LineWidth',5,'Marker','o','MarkerFaceColor','r','MarkerSize',15);
plot(nanmedian(x2(:,:),1)/nansum(nanmedian(x2(:,:),1))*100,'Color','b','LineWidth',5,'Marker','o','MarkerFaceColor','b','MarkerSize',15);


l = legend(["ipsi","contra"]);
set(l,'Color','none');

xlim([1,30]);
ylim([0,25]);
xtick = get(gca,'XTick');
xticks(1:2:max(xtick));
xticklabels(1:2:max(xtick));
xtickangle(45);
xlabel('N');
ylabel('P(%)');

set(gca,'Fontsize',25,'Color','none');
set(gcf,'Position',[0,0,1000,800]);



%% Fig. 2F
for i = 1:length(PFCregions)
    figure
    hold on
    p1 = plot(nanmean(x1(i,:),1)/nansum(nanmean(x1(i,:),1)),'Color','r','LineWidth',5,'Marker','o','MarkerFaceColor','r','MarkerSize',15);
    p2 = plot(nanmean(x2(i,:),1)/nansum(nanmean(x2(i,:),1)),'Color','b','LineWidth',5,'Marker','o','MarkerFaceColor','b','MarkerSize',15);


    regionname = PFCregions(i);
    title(regionname);
    l = legend([p1,p2],["ipsi","contra"]);
    set(l,'Color','none');
    set(gca,'Fontsize',15,'Color','none');
    xlim([1,30]);
    xtick = get(gca,'XTick');
    xticks(1:2:max(xtick));
    xticklabels(1:2:max(xtick));
    xtickangle(45);
%     xlabel('N');
%     ylabel('P(%)');
    set(gcf,'Position',[0,0,1000,800]);


    figname = [char(regionname),'_bilateral_divergence'];

end

%%





















