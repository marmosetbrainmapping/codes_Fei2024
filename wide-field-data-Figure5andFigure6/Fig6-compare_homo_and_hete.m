clear;
clc;
close all;
% data1_path = '/wide-field-data/structure_connectivity/R2L_SC_strength.mat';
% data1 = load(data1_path).connectivity;
% data1(data1<10.^-1.5) = 0;

data2_path = '/wide-field-data/contraFC_mean.mat';
data2 = load(data2_path).contraFC_mean;
SC_path = '/wide-field-data/structure_connectivity/full_SC_mask_wf.mat';
SC_mask = load(SC_path).full_SC_mask_wf;   
contra_mask = SC_mask(1:31,32:62);
data2 = data2 .* contra_mask;

% SC_path ='/wide-field-data/structure_connectivity/full_SC_mask_wf.mat';
% SC_mask =  load(SC_path).full_SC_mask_wf;  
% [hete_idx, hete_vec] = get_hete_FC(SC_mask);
% 
% saliences = load('./new_analysis/condition_varibility/saliences_conn_type3.mat').saliences_mean;
% saliences = 1 ./ saliences;
% contralateral = zeros(31,31);
% contralateral(hete_idx{2,1}) = saliences(1454+726/2+1: 1454+726);
% 
% homo = saliences(1454+726+1:2242);
% for i=1:31
%     contralateral(i,i) = homo(i);
% end
% data3 = contralateral .* SC_mask(1:31,32:62);



% region names:31
fileID = fopen('/wide-field-data/preprocess/ccfv3/mask_label_name_correct.txt');
file = textscan(fileID, '%s');
fclose(fileID);
regionNames=file{1,1};
%region names:43
fileID1 = fopen('/wide-field-data/structure_connectivity/cortex_region_names.txt');
file1 = textscan(fileID1, '%s');
fclose(fileID1);
regionNames1 = file1{1,1};

%module region names:
fileID2 = fopen('/wide-field-data/structure_connectivity/cortex_name_module.txt');
cortex_module_name = textscan(fileID2,'%s');
fclose(fileID2);
cortex_module_name = cortex_module_name{1,1};

order_43 = {};   % the relation 
for i =1:43
    n = cortex_module_name{i,1};
    for j =1:43
        if strcmp(regionNames1{j,1},n)
            order_43{end+1} = j;
        end
    end
end
order_43 = cell2mat(order_43');

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




region_names = regionNames(order_31);
data = data2(order_31, order_31);

homo = diag(data);

mask_x = {}; % hete > homo
mask_y = {};
region_y = {};

data_points_x = {};
data_points_y = {};
for i=1:31
    row = data(i,:);
    row(i) = 0;
    idx = find(row>0);
    hete = row(idx);
    for j=1:length(hete)
        data_points_x{end+1}=i;
        data_points_y{end+1}=hete(j);
        if hete(j)>homo(i)
            mask_x{end+1} = i;
            mask_y{end+1} = hete(j);
            region_y{end+1} = region_names(find(row==hete(j)));
        end
    end
end
data_points_x = cell2mat(data_points_x);
data_points_y = cell2mat(data_points_y);
mask_x = cell2mat(mask_x);
mask_y = cell2mat(mask_y);

x = 1:1:31;
bar(x, homo);
% scatter(x, homo, 50, 'Marker', '*')
ylim([0.3 1]);
x = categorical(region_names);
x = reordercats(x,region_names);
xticks([]);
xticks(1:31);
xticklabels(region_names);
xtickangle(90);
hold on;
scatter(data_points_x, data_points_y,40, 'filled');
for i=1:length(mask_x)
    text(mask_x(i), mask_y(i), region_y{i}, 'HorizontalAlignment', 'center','VerticalAlignment','middle')
end


function [hete_idx, hete_FC_vec] = get_hete_FC(matrix)
matrix1 = matrix(1:31,32:62);
matrix2 = matrix(32:62,1:31);
lower1 = tril(matrix1,-1);
upper1 = triu(matrix1,1);
lower2 = tril(matrix2,-1);
upper2 = triu(matrix2,1);
nonzero_idx1 = find(lower1);
nonzero_idx2 = find(upper1);
nonzero_idx3 = find(lower2);
nonzero_idx4 = find(upper2);
% nonzero_idx = [nonzero_idx1;nonzero_idx2];
hete_FC_vec1 = lower1(nonzero_idx1);
hete_FC_vec2 = upper1(nonzero_idx2);
hete_FC_vec3 = lower2(nonzero_idx3);
hete_FC_vec4 = upper2(nonzero_idx4);
hete_FC_vec = [hete_FC_vec1; hete_FC_vec2;hete_FC_vec3;hete_FC_vec4];
hete_idx = {};
hete_idx{1,1} = [nonzero_idx1; nonzero_idx2];
hete_idx{2,1} = [nonzero_idx3; nonzero_idx4];
end

