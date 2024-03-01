clear;
clc;
close all;
% SC_path ='/wide-field-data/structure_connectivity/full_SC_type5_mask_wf.mat';
% sc_mask =  load(SC_path).full_SC_type5_mask_wf;
sc_path = '/wide-field-data/structure_connectivity/full_SC_mask_wf.mat';
sc_mask = load(sc_path).full_SC_mask_wf;

%% get FC strength of every connection in ipsi and contra
ori_path = '/wide-field-data/valid_data/';
% ori_path = '/wide-field-data/new_analysis/distinguish_different_state_result/results/';
% ori_path = '/wide-field-data/filter/data2/';
fileID = fopen('/wide-field-data/new_analysis/distinguish_different_state_result/present_sessions.txt');
sessions = textscan(fileID, '%s');
fclose(fileID);
sessions_name = sessions{1,1};
sample_num = length(sessions_name);

contraFC_sum = zeros(31,31);
ipsiFC_sum = zeros(31,31);

for i=1:sample_num
    sessionName = char(sessions_name(i));
    data_path = fullfile(ori_path, sessionName);
    data_path = char(data_path);
    file_path = dir(fullfile(data_path,'full_FC.mat'));
    if isempty(file_path)
        fprintf('Wrong file path!\n');
        return;
    end 
    file_path = fullfile(file_path(1).folder, file_path(1).name);
    full_FC=load(file_path).full_FC;
    full_FC=full_FC.*sc_mask;
    ipsi_FC=full_FC(1:31,1:31);
    contra_FC=full_FC(1:31,32:62);
    ipsiFC_sum=ipsiFC_sum+ipsi_FC;
    contraFC_sum=contraFC_sum+contra_FC;
end
contraFC_mean = contraFC_sum / sample_num;
ipsiFC_mean = ipsiFC_sum / sample_num;

%% get 3 type relation: ipsi exists and hete unexists , ipsi FC > hete FC, ipsi FC < hete FC
fileID = fopen('/wide-field-data/preprocess/ccfv3/mask_label_name_correct.txt');
file = textscan(fileID, '%s');
fclose(fileID);
regionNames=file{1,1};

%module region names:
fileID2 = fopen('/wide-field-data/structure_connectivity/cortex_name_module.txt');
cortex_module_name = textscan(fileID2,'%s');
fclose(fileID2);
cortex_module_name = cortex_module_name{1,1};

order_31 = {};   % the relation 
for i =1:43
    n = cortex_module_name{i,1};
    for j =1:31
        if strcmp(regionNames{j,1},n)
            order_31{end+1} = j;
        end
    end
end
order_31 = cell2mat(order_31');

ipsiFC_mean = ipsiFC_mean(order_31, order_31);
contraFC_mean = contraFC_mean(order_31, order_31);

ipsi = sc_mask(1:31,1:31);  %ipsi  structure connectivity mask
contra = sc_mask(1:31,32:62); %contra structure connectivity mask
ipsi = ipsi(order_31, order_31);
contra = contra(order_31, order_31);
source_table = {};
target_table = {};
value_table = {};

for i=1:31
    for j=1:31
        if i == j
            continue;
        end
      source_table{end+1} = i ;  
      target_table{end+1} = j ;
      if contra(i,j) ~= 0 & ipsi(i,j) ~= 0
%           value_table{end+1} = 1;   %ipsi and contra  all exist
           if ipsiFC_mean(i,j) > contraFC_mean(i,j)
              value_table{end+1}=2;  % ipsi > hete
           elseif ipsiFC_mean(i,j) < contraFC_mean(i,j)
              value_table{end+1}=3;  % ipsi <= hete
          end
      elseif contra(i,j) == 0 & ipsi(i,j) ~= 0
          value_table{end+1} = -1;  % ipsi exists but contra unexists
      elseif contra(i,j) ~=0 & ipsi(i,j) == 0
          value_table{end+1} = -2; % ipsi unexists but contra exists [this sitution may be unexisting]
      else
          value_table{end+1} = 0;
      end
          
        
    end
end

source_table = source_table';
target_table = target_table';
value_table = value_table';

value = cell2mat(value_table);
source = cell2mat(source_table);
target = cell2mat(target_table);
idx_1 = find(value==-1);
source_1 = source(idx_1);
target_1 = target(idx_1);

idx_2 = find(value==2);
source_2 = source(idx_2);
target_2 = target(idx_2);

idx_3 = find(value==3);
source_3 = source(idx_3);
target_3 = target(idx_3);

figure;
scatter(source_1, target_1, [],[0 0.376 0.701]);
xticks(1:31);
xticklabels(regionNames(order_31));
xtickangle(90);
yticks(1:31);
yticklabels(regionNames(order_31));
hold on;
scatter(source_2,target_2,[],[0.569 0.824 0.753],'filled');
hold on;
scatter(source_3,target_3,[],[0.7373 0.2431 0.012],'*');
legend('only intra','intra>hete','intra<hete');
