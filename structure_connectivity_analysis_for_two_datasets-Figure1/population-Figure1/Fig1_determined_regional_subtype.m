clear;
clc;
close all;
addpath('/mnt/Data16Tb/Data/feiyao/wield-field-data/preprocess/code/colormap');
data_path = '/mnt/Data16Tb/Data/feiyao/wield-field-data/structure_connectivity/surface_parcellation.mat';
data1_path = '/mnt/Data16Tb/Data/feiyao/wield-field-data/structure_connectivity/volume_SC/new_sparse/R2L_SC_strength.mat';
data2_path = '/mnt/Data16Tb/Data/feiyao/wield-field-data/structure_connectivity/volume_SC/new_sparse/R2R_SC_strength.mat';
region_name_path = '/mnt/Data16Tb/Data/feiyao/wield-field-data/structure_connectivity/volume_SC/cortex_region_names.txt';
fileID = fopen(region_name_path);
region_name = textscan(fileID,'%s');
fclose(fileID);
region_name = region_name{1,1};
region_num = length(region_name);

fileID2 = fopen('cortex_name_module.txt');
cortex_module_name = textscan(fileID2,'%s');
fclose(fileID2);
cortex_module_name = cortex_module_name{1,1};
 
order = {};   % the relation 
for i =1:43
    n = cortex_module_name{i,1};
    for j =1:43
        if strcmp(region_name{j,1},n)
            order{end+1} = j;
        end
    end
end
order = cell2mat(order');
r = region_name(order);
% order = flip(order);
% cortex_module_name = flip(cortex_module_name);

R2L_SC = load(data1_path).connectivity;
R2R_SC = load(data2_path).connectivity;
thresh = 10 .^ -1.5;
R2L_SC(R2L_SC<thresh)=0;
R2R_SC(R2R_SC<thresh)=0;
% full_SC = load('full_SC_type5.mat').full_SC_type5;  %  cell threshold
% R2L_SC = full_SC(44:86, 1:43);
% R2R_SC = full_SC(44:86, 44:86);

% calculate specific ipsilateral region, specific contralateral region and
% common region
b = {};  % common region
si = {}; % specific ipsilateral region
sc = {}; % specific contralateral region
R2R_SC = R2R_SC(order, order);
R2L_SC = R2L_SC(order, order);
for i=1:43
% R2R_SC(i, i) = 0;
iii = R2R_SC(i,:);
ccc = R2L_SC(i,:);
iii_non = find(iii);
ccc_non = find(ccc);
b{end+1} = length(intersect(iii_non, ccc_non));
si{end+1} = length(setdiff(iii_non, ccc_non)); 
sc{end+1} = length(setdiff(ccc_non,iii_non));
end
symmetry = horzcat(cell2mat(si)', cell2mat(b)', cell2mat(sc)');  % primary threshold:37 IB, 6 IBC
                                                                                                          % cell threshold: 41 IB, 2 IBC

 
b_nonzero_idx = find(symmetry(:, 2));
c_ = symmetry(:,3);
c_ = c_(b_nonzero_idx);
figure;
sizes = [length(find(c_==0))/length(b_nonzero_idx), length(find(c_~=0))/length(b_nonzero_idx)];
labels = {'IB', 'IBC'};
pie(sizes, labels);
saveas(gcf, '/mnt/Data16Tb/Data/feiyao/wield-field-data/structure_connectivity/volume_SC/new_sparse/regional_types/WT_percentage.pdf');
