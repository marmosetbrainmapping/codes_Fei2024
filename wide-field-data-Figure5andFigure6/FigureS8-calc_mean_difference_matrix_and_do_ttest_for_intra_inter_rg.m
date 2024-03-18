clear;
clc;
close all;
addpath('/mnt/Data16Tb/Data/feiyao/wield-field-data/preprocess/code/colormap/')

region_name_path = '/mnt/Data16Tb/Data/feiyao/wield-field-data/preprocess/ccfv3/mask_label_name_correct.txt';
fileID = fopen(region_name_path);
region_name = textscan(fileID,'%s');
fclose(fileID);
region_name = region_name{1,1};
region_num = length(region_name);

fileID2 = fopen('/mnt/Data16Tb/Data/feiyao/wield-field-data/preprocess/ccfv3/mask_label_name_correct_reorder_.txt');
cortex_module_name = textscan(fileID2,'%s');
fclose(fileID2);
cortex_module_name = cortex_module_name{1,1};
 
order = {};   % the relation 
for i =1:31
    n = cortex_module_name{i,1};
    for j =1:31
        if strcmp(region_name{j,1},n)
            order{end+1} = j;
        end
    end
end
order = cell2mat(order');

% ori_path = '/mnt/Data16Tb/Data/feiyao/wield-field-data/valid_data/';
ori_path = '/mnt/Data16Tb/Data/feiyao/wield-field-data/new_analysis/distinguish_different_state_result/results/';
fileID = fopen('/mnt/Data16Tb/Data/feiyao/wield-field-data/new_analysis/distinguish_different_state_result/present_sessions.txt');
sessions = textscan(fileID, '%s');
fclose(fileID);
sessions_name = sessions{1,1};
sample_num = length(sessions_name);


intra_inter_all = zeros(31,31);
statistic_intra_inter = zeros(31*31, sample_num);
for i=1:sample_num
    sessionName = char(sessions_name(i));
    data_path = fullfile(ori_path, sessionName);
    data_path = char(data_path);
    file_path = fullfile(data_path,'rem_intra-inter_FC_rg.mat');  %2:L   1:R
    difference = load(file_path).difference;
    intra_inter_all = intra_inter_all + difference;
    difference_vector = reshape(difference, [], 1);
    statistic_intra_inter(:,i) = difference_vector;
end

p_value_difference = zeros(31*31, 1);
pp_value_difference = zeros(31*31, 1);
stats_difference = zeros(31*31, 1);
for j=1:961
    vector = statistic_intra_inter(j, :);
    % 执行单样本 t 检验
    [h, p, ci, stats] = ttest(vector', 0);
    
    p_value_difference(j) = h;
    stats_difference(j) = stats.tstat;
    pp_value_difference(j) = p;
end

p_difference = reshape(p_value_difference, 31, 31);
p_difference(isnan(p_difference)) = 0;
stats_difference = reshape(stats_difference, 31, 31);
stats_difference = stats_difference(order, order);
pp_value_difference = reshape(pp_value_difference, 31, 31);
pp_value_difference = pp_value_difference(order, order);

difference_mean = intra_inter_all / sample_num;

difference_mean = difference_mean .* p_difference;

difference_mean = difference_mean(order,order);

ratio = length(find(difference_mean >0)) / (31*31);
ratio2 = length(find(difference_mean >0)) / (31*31);
ratio1 = (length(find(difference_mean >0))-31) / (31*30);

% RdBu = flipud(RdBu(127));
% map = RdBu;
% map(64, :) = [137/255 137/255 137/255];
% figure;
% set(gcf, 'Position', [0, 0, 600, 520]);
% imagesc(difference_mean);
% xticks(1:31);
% yticks(1:31);
% xticklabels(cortex_module_name);
% xtickangle(90);
% yticklabels(cortex_module_name);
% caxis([-max(difference_mean(:)), max(difference_mean(:))]);
% colormap(map);
% colorbar;
% saveas(gcf, '/mnt/Data16Tb/Data/feiyao/wield-field-data/difference_results/new/rem_intra_inter_rg.pdf');
% 
% origin =  intra_inter_all / sample_num;
% origin = origin(order,order);
% figure;
% set(gcf, 'Position', [0, 0, 600, 520]);
% imagesc(origin);
% xticks(1:31);
% yticks(1:31);
% xticklabels(cortex_module_name);
% xtickangle(90);
% yticklabels(cortex_module_name);
% caxis([-max(origin(:)), max(origin(:))]);
% colormap(RdBu);
% colorbar;
% saveas(gcf, '/mnt/Data16Tb/Data/feiyao/wield-field-data/difference_results/new/rem_intra_inter_rg_origin.pdf');
