clear;
clc;
close all;
ori_path = '/wide-field-data/new_analysis/distinguish_different_state_result/results/';
fileID = fopen('/wide-field-data/new_analysis/distinguish_different_state_result/present_sessions.txt');
sessions = textscan(fileID, '%s');
fclose(fileID);
sessions_name = sessions{1,1};
sample_num = length(sessions_name);
% SC_path = '/wide-field-data/structure_connectivity/full_SC_mask_wf.mat';
% SC_mask = load(SC_path).full_SC_mask_wf;
SC_path ='/wide-field-data/structure_connectivity/full_SC_type5_mask_wf.mat';
SC_mask =  load(SC_path).full_SC_type5_mask_wf;

homoFC_all_subs = zeros(sample_num,1);
heteFC_all_subs = zeros(sample_num,1);
ipsiFC_all_subs = zeros(sample_num,1);

for i=1:sample_num
    sessionName = char(sessions_name(i));
    data_path = fullfile(ori_path, sessionName);
    data_path = char(data_path);
    file_path = dir(fullfile(data_path,'wake_full_FC.mat'));
    if isempty(file_path)
        fprintf('Wrong file path!\n');
        return;
    end 
    file_path = fullfile(file_path(1).folder, file_path(1).name);
   % load(file_path, 'homo_FC_vector', 'hete_FC_vector', 'ipsi_FC_vector');
    full_FC = load(file_path).full_FC;
    full_FC = full_FC.* SC_mask;

  
     
      homo_FC_vector = get_ho_FC(full_FC);
      hete_FC_vector = get_he_FC(full_FC);
      ipsi_FC_vector = get_ipsi_FC(full_FC);

    FC_all(i,:) = [ipsi_FC_vector; hete_FC_vector; homo_FC_vector];
    homo_mean = mean(homo_FC_vector);
    homoFC_all_subs(i) = homo_mean;
    hete_mean = mean(hete_FC_vector);
    heteFC_all_subs(i) = hete_mean;
    ipsi_mean = mean(ipsi_FC_vector);
    ipsiFC_all_subs(i) = ipsi_mean; 
end
FC_mean = mean(FC_all, 1)';

FC = [ipsiFC_all_subs, heteFC_all_subs, homoFC_all_subs];
FC_mean = mean(FC,1);
FC_sem = std(FC,1) / sqrt(sample_num);
figure;
X = categorical({'ipsi','hete','homo'});
X = reordercats(X,{'ipsi','hete','homo'});
b = bar(X,FC_mean, 0.5);
ylim([0,1]);
b.FaceColor = 'flat';
b.CData(1,:) = [0.902 0.489 0.482];
b.CData(2,:) = [0.5608 0.7647 0.1216];
b.CData(3,:) = [0.4157 0.6314 1];
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',7);
ylabel('FC')
xlabel('Mean+-SEM')
hold on;
plot([X(1), X(1)], [FC_mean(1) - FC_sem(1), FC_mean(1) + FC_sem(1)], 'color',[0 0 0],'Linestyle', '-');
hold on;
plot([X(2), X(2)], [FC_mean(2) - FC_sem(2), FC_mean(2) + FC_sem(2)], 'color',[0 0 0],'Linestyle', '-');
hold on;
plot([X(3), X(3)], [FC_mean(3) - FC_sem(3), FC_mean(3) + FC_sem(3)], 'color',[0 0 0],'Linestyle', '-');
% hold on;
% scatter(repmat(X(1), sample_num, 1), FC(:,1),[],[0.4471 0.4441 0.4441]);
% hold on;
% scatter(repmat(X(2), sample_num, 1), FC(:,2),[],[0.4471 0.4441 0.4441]);
% hold on;
% scatter(repmat(X(3), sample_num, 1), FC(:,3),[],[0.4471 0.4441 0.4441]);



% normal distribution test and variance 
[~,col]=size(FC);
h_test = zeros(col,1);
for i=1:col
    h=lillietest(FC(:,i),'alpha',0.001);  % return 0, normal distribution
    h_test(i)=h;
end
if length(find(h_test==0)) >= 1
   if vartestn(FC) > 0.001
       [p,tb1,stats] = anova1(FC);
   else
      fprintf('variances are not same\n');
      [p, tb1, stats] = kruskalwallis(FC);
   end
else
    fprintf('distributions are not normal\n');
   [p, tb1, stats] = kruskalwallis(FC);
end

c = multcompare(stats);
tb2 = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

function ho_FC_vec = get_ho_FC(matrix)
ho_FC1 = diag(matrix(1:31, 32:62));
ho_FC2 = diag(matrix(32:62, 1:31));
nonzero_idx1 = find(ho_FC1);  
nonzero_idx2 = find(ho_FC2);  
ho_FC_vec1 = ho_FC1(nonzero_idx1);
ho_FC_vec2 = ho_FC2(nonzero_idx2);
ho_FC_vec = [ho_FC_vec1;ho_FC_vec2];
end

function ipsi_FC_vec = get_ipsi_FC(matrix)
matrix1 = matrix(1:31,1:31);
matrix2 = matrix(32:62,32:62);
lower1 = tril(matrix1,-1);
upper1 = triu(matrix1,1);
lower2 = tril(matrix2,-1);
upper2 = triu(matrix2,1);
nonzero_idx1 = find(lower1);
nonzero_idx2 = find(upper1);
nonzero_idx3 = find(lower2);
nonzero_idx4 = find(upper2);
ipsi_FC_vec1 = lower1(nonzero_idx1);
ipsi_FC_vec2 = upper1(nonzero_idx2);
ipsi_FC_vec3 = lower2(nonzero_idx3);
ipsi_FC_vec4 = upper2(nonzero_idx4);
ipsi_FC_vec = [ipsi_FC_vec1; ipsi_FC_vec2;ipsi_FC_vec3;ipsi_FC_vec4];
end

function hete_FC_vec = get_he_FC(matrix)
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
end



