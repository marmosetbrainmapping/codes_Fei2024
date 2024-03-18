load('ipsi_neuron_to_region_count.mat');
ipsi_count = counts;
load('contra_neuron_to_region_count.mat');
contra_count = counts;

% ipsi_count & contra_count both have size 3704*43.
% Each row represents a neuron, each column represents a cortical region.
% The value at i row j colum represents the number of projection terminals of the projections from neuron i to region j.
%%
load('ITcortexneuronname.mat');
load('cortexregionnames.mat');

PFCregions = unique(ITcortexneuronname);
Cortexregions = (cortexregionnames);

PFCinds = zeros(length(PFCregions),1);
for i = 1:length(PFCregions)
    PFCregion = PFCregions(i);
    ind = find(Cortexregions==PFCregion);
    PFCinds(i) = ind;
end

ipsi_count_pfc = ipsi_count(:,PFCinds);
contra_count_pfc = contra_count(:,PFCinds);
neuroninds = cell(length(PFCinds),1);
for i = 1:length(PFCinds)
    PFCregion = PFCregions(i);
    inds = find(ITcortexneuronname==PFCregion);
    neuroninds{i,1} = inds;
end
% ipsi_count_pfc & contra_count_pfc both have size 3704*11, representing the projections to 11 PFC regions.


ipsi_count_neurons = cell(length(PFCinds),1);
contra_count_neurons = cell(length(PFCinds),1);
for i = 1:length(PFCinds)
    inds = neuroninds{i};
    ipsi_count_neurons{i,1} = ipsi_count(inds,:);
    contra_count_neurons{i,1} = contra_count(inds,:);
end
% ipsi_count_neurons & contra_count_neurons are both 11x1 cells. Each row
% element of the cell represents the projections of neurons within a
% specific PFC regions to all 43 cortical regions.


ipsi_proj_neurons = cell(length(PFCregions),1);
contra_proj_neurons = cell(length(PFCregions),1);
for i = 1:length(PFCregions)

    i_projs = ipsi_count_neurons{i};
    ipsi_proj_neurons{i} = double(i_projs>0);

    i_projs = contra_count_neurons{i};
    contra_proj_neurons{i} = double(i_projs>0);

end
% ipsi_proj_neurons & contra_proj_neurons are both 11x1 cells. Each row
% element of the cell is the binarized result of the corresponding element
% of ipsi_count_neurons & contra_count_neurons.
%%
data_heterogeneity = zeros(length(PFCregions),length(Cortexregions));

nprojections = zeros(length(PFCregions),length(Cortexregions));

for i = 1:size(data_heterogeneity,1)

    ipsi_i_projs = ipsi_proj_neurons{i};
    contra_i_projs = contra_proj_neurons{i};

    ipsi_i_terminals = ipsi_count_neurons{i};
    contra_i_terminals = contra_count_neurons{i};
    
    n_neuronipsi = size(ipsi_i_projs,1);
    n_neuroncontra = size(contra_i_projs,1);

    for j = 1:size(data_heterogeneity,2)
    
        ipsi_i2j_projs = ipsi_i_projs(:,j);
        contra_i2j_projs = contra_i_projs(:,j);

        ipsi_inds = find(ipsi_i2j_projs==1);
        contra_inds = find(contra_i2j_projs==1);
        common_inds = intersect(ipsi_inds,contra_inds);
        union_inds = union(ipsi_inds,contra_inds);
        unique_inds = setdiff(union_inds,common_inds);
        ipsi_unique_inds = setdiff(ipsi_inds,common_inds);
        contra_unique_inds = setdiff(contra_inds,common_inds);

        p1 = length(ipsi_inds)/n_neuronipsi;
        p2 = length(contra_inds)/n_neuroncontra;

        overlap = length(common_inds)/(min(length(ipsi_inds),length(contra_inds)));
        data_heterogeneity(i,j) = 1-overlap;
        nprojections(i,j) = length(ipsi_inds) + length(contra_inds);

    end

end

% data_heterogeneity has the size 11x43. Each row represents an upstream
% region of projections (among the 11 PFC regions). Each column represents
% a downstream region (among the 43 cortical regions). The value at i row j
% column represents the heterogeneity of projections from the PFC region i
% to a pair of homotopic downstream regions j.
%%
%%
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
%%
cmap = colormap(cbrewer2('RdBu', 100));
cmap = flipud(cmap);
%%
figure
hold on

x = data_heterogeneity(PFC_dataname2figname,cortex_dataname2figname);
nx = length(cortex_dataname2figname);
ny = length(PFC_dataname2figname);
cdata = zeros(size(x,1),size(x,2),3);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        c = x(i,j);
        if ~isnan(c)
            ind = round(c*size(cmap,1));
            if ind == 0
                ind = 1;
            elseif ind > size(cmap,1)
                ind = size(cmap,1);
            end
            cdata(i,j,:) = cmap(ind,:);
        else
            cdata(i,j,:) = [1,1,1] * 0.5;
        end
    end
end
imagesc('CData',flipud(cdata));
axis equal
axis tight

xticks(1:length(fig_cortexregions));
xticklabels(fig_cortexregions);
yticks(1:length(fig_PFCregions));
yticklabels(fig_PFCregions);

set(gca,'Fontsize',15);
xtickangle(45);
set(gca,'xtick', linspace(0.5,nx+0.5,nx+1),'ytick', linspace(0.5,ny+0.5,ny+1));
set(gca,'xgrid','on','ygrid','on','GridColor',[1,1,1]*0.9,'gridlinestyle','-');
set(gcf,'Position',[1036 850 1492 480]);
