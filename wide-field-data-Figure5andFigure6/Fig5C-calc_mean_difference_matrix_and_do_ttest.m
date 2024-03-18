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

fileID2 = fopen('/mnt/Data16Tb/Data/feiyao/wield-field-data/preprocess/ccfv3/mask_label_name_correct_reorder.txt');
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

ori_path = '/mnt/Data16Tb/Data/feiyao/wield-field-data/new_analysis/distinguish_different_state_result/results/';
fileID = fopen('/mnt/Data16Tb/Data/feiyao/wield-field-data/new_analysis/distinguish_different_state_result/present_sessions.txt');
sessions = textscan(fileID, '%s');
fclose(fileID);
sessions_name = sessions{1,1};
sample_num = length(sessions_name);

intra_all = zeros(31,31);
inter_all = zeros(31,31);
statistic_intra = zeros(31*31, sample_num);
statistic_inter = zeros(31*31, sample_num);
for i=1:sample_num
    sessionName = char(sessions_name(i));
    data_path = fullfile(ori_path, sessionName);
    data_path = char(data_path);
    file_path = fullfile(data_path,'rem-nrem_FC.mat');
    difference = load(file_path).difference;
    difference_R2R = difference(32:62, 32:62);
    difference_R2L = difference(32:62,1:31);
    intra_all = intra_all + difference_R2R;
    inter_all =  inter_all + difference_R2L;
    difference_R2R_vector = reshape(difference_R2R, [], 1);
    difference_R2L_vector = reshape(difference_R2L, [], 1);
    statistic_intra(:,i) = difference_R2R_vector;
    statistic_inter(:,i) = difference_R2L_vector;
    
end

p_value_intra = zeros(31*31, 1);
p_value_inter = zeros(31*31, 1);
for j=1:961
    intra_vector = statistic_intra(j, :);
    inter_vector = statistic_inter(j, :);
    % 执行单样本 t 检验
    [h, p, ci, stats] = ttest(intra_vector', 0);
    [h1, p1, ci1, stats1] = ttest(inter_vector', 0);
    
    p_value_intra(j) = h;
    p_value_inter(j) = h1;
end

p_intra = reshape(p_value_intra, 31, 31);
p_inter = reshape(p_value_inter, 31, 31);
p_intra(isnan(p_intra)) = 0;

intra_mean = intra_all / sample_num;
inter_mean = inter_all/ sample_num;

intra_mean = intra_mean .* p_intra;
inter_mean = inter_mean .* p_inter;

intra_mean = intra_mean(order,order);
inter_mean = inter_mean(order,order);

% map = coolwarm(127);
% map(64, :) = [137/255 137/255 137/255];
% figure;
% set(gcf, 'Position', [0, 0, 600, 520]);
% imagesc(intra_mean);
% xticks(1:31);
% yticks(1:31);
% xticklabels(cortex_module_name);
% xtickangle(90);
% yticklabels(cortex_module_name);
% caxis([-max(intra_mean(:)), max(intra_mean(:))])
% colormap(map);
% colorbar;
% saveas(gcf, '/mnt/Data16Tb/Data/feiyao/wield-field-data/difference_results/rem-nrem_intra.pdf');
% 
% figure;
% set(gcf, 'Position', [0, 0, 600, 520]);
% imagesc(inter_mean);
% xticks(1:31);
% yticks(1:31);
% xticklabels(cortex_module_name);
% xtickangle(90);
% yticklabels(cortex_module_name);
% caxis([-max(inter_mean(:)), max(inter_mean(:))])
% colormap(map);
% colorbar;
% saveas(gcf, '/mnt/Data16Tb/Data/feiyao/wield-field-data/difference_results/rem-nrem_inter.pdf');
% 
% 

    
   