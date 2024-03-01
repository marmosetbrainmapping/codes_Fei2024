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
%
ipsi_proj_neurons = cell(length(PFCregions),1);
contra_proj_neurons = cell(length(PFCregions),1);
%
for i = 1:length(PFCregions)

    i_projs = ipsi_count_neurons{i};
    ipsi_proj_neurons{i} = double(i_projs>0);

    i_projs = contra_count_neurons{i};
    contra_proj_neurons{i} = double(i_projs>0);

end

fig_PFCregions = ["FRP","ACAd","ACAv","PL","ILA","ORBl","ORBm","ORBvl","AId","AIv","MOs"];
fig_cortexregions = ["FRP","ACAd","ACAv","PL","ILA","ORBl","ORBm","ORBvl","AId","AIv","MOs","AIp","GU","VISC","TEa","PERI","ECT","SSs","SSp-bfd","SSp-tr","SSp-ll","SSp-ul","SSp-un","SSp-n","SSp-m","MOp","VISal","VISl","VISp","VISpl","VISli","VISpor","VISrl","VISa","VISam","VISpm","RSPagl","RSPd","RSPv","AUDd","AUDp","AUDpo","AUDv"];
%
PFC_dataname2figname = zeros(length(PFCregions),1);
cortex_dataname2figname = zeros(length(cortexregionnames),1);

for i = 1:length(PFCregions)
    figname = fig_PFCregions(i);
    dataind = find(PFCregions == figname);
    PFC_dataname2figname(i) = dataind;
end

for i = 1:length(cortexregionnames)
    figname = fig_cortexregions(i);
    dataind = find(cortexregionnames == figname);
    cortex_dataname2figname(i) = dataind;
end

load('ccfv3bilateralFC.mat');
load('ccfv3cortexregionnames.mat');
[intersect_regions,inds,indsccfv3] = intersect(cortexregionnames,ccfv3cortexregionnames);
regionscorr = ccfv3bilateralFC(indsccfv3);
%%
heterogeneity_downstream = zeros(1,length(cortexregionnames));
heterogeneity_upstream = zeros(1,length(PFCregions));

for i = 1:length(heterogeneity_downstream)

    ipsi_projs = [];
    contra_projs = [];
    for j = 1:length(PFCregions)
        if PFCregions(j) == cortexregionnames(i)
            continue
        else
            ipsi_projs = [ipsi_projs;ipsi_count_neurons{j}];
            contra_projs = [contra_projs;contra_count_neurons{j}];
        end
    end
    ipsi_inds = find(ipsi_projs(:,i)~=0);
    contra_inds = find(contra_projs(:,i)~=0);
    common_inds = intersect(ipsi_inds,contra_inds);
    overlap = length(common_inds)/min(length(ipsi_inds),length(contra_inds));
    heterogeneity_downstream(i) = 1-overlap;

end

for i = 1:length(heterogeneity_upstream)

    ipsi_i_projs = ipsi_proj_neurons{i};
    contra_i_projs = contra_proj_neurons{i};
    ipsi_inds = find(sum(ipsi_i_projs,2)~=0);
    contra_inds = find(sum(contra_i_projs,2)~=0);
    common_inds = intersect(ipsi_inds,contra_inds);
    overlap = length(common_inds)/min(length(ipsi_inds),length(contra_inds));
    heterogeneity_upstream(i) = 1-overlap;

end
%%
figure
plot(heterogeneity_downstream(cortex_dataname2figname),'LineWidth',3,'Color','black','LineStyle','-');
xticks(1:length(heterogeneity_downstream));
xticklabels(fig_cortexregions);
ylabel('heterogeneity');
set(gca,'Color','none');
set(gcf,'Position',[1036 850 1492 250]);
set(gca,'Fontsize',15);
xtickangle(45);



