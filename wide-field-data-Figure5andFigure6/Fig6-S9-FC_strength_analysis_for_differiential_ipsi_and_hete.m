clear;
clc;
close all;
ori_path = '/wide-field-data/valid_data/';
% ori_path = '/wide-field-data/new_analysis/distinguish_different_state_result/results/';
% SC_path = '/wide-field-data/structure_connectivity/full_SC_mask_wf.mat';
% SC_mask = load(SC_path).full_SC_mask_wf;
SC_path ='/wide-field-data/structure_connectivity/full_SC_type5_mask_wf.mat';
SC_mask =  load(SC_path).full_SC_type5_mask_wf;
fileID = fopen('/wide-field-data/new_analysis/distinguish_different_state_result/present_sessions.txt');
sessions = textscan(fileID, '%s');
fclose(fileID);
sessions_name = sessions{1,1};
sample_num = length(sessions_name);

heteFC_all_subs = zeros(sample_num,1);
ipsiFC_all_subs = zeros(sample_num,1);
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
    full_FC = load(file_path).full_FC;
    full_FC = full_FC .* SC_mask;
    [ipsi_FC_vector, hete_FC_vector] = maintain_specific_connections(full_FC);
    hete_mean = mean(hete_FC_vector);
    heteFC_all_subs(i) = hete_mean;
    ipsi_mean = mean(ipsi_FC_vector);
    ipsiFC_all_subs(i) = ipsi_mean; 
    
end

FC = [ipsiFC_all_subs, heteFC_all_subs];
FC_mean = mean(FC,1);
FC_sem = std(FC,1) / sqrt(sample_num-1);
figure;
X = categorical({'ipsi','hete'});
X = reordercats(X,{'ipsi','hete'});
b = bar(X,FC_mean, 0.5);
ylim([0 1]);
b.FaceColor = 'flat';
b.CData(1,:) = [0.035 0.576 0.588];
b.CData(2,:) = [0.682 0.125 0.071];
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',7);
ylabel('FC Strength')
xlabel('Mean+-SEM')
hold on;
plot([X(1), X(1)], [FC_mean(1) - FC_sem(1), FC_mean(1) + FC_sem(1)], 'color',[0.2, 0.2, 0.2],'Linestyle', '-');
hold on;
plot([X(2), X(2)], [FC_mean(2) - FC_sem(2), FC_mean(2) + FC_sem(2)], 'color',[0.2, 0.2, 0.2],'Linestyle', '-');




% normal distribution test and variance 
[~,col]=size(FC);
h_test = zeros(col,1);
for i=1:col
    h=lillietest(FC(:,i));
    h_test(i)=h;
end
if length(find(h_test==0)) >= 1
   if vartestn(FC) > 0.05
       [p,tb1,stats] = anova1(FC);
   else
      fprintf('variances are not same\n');
      [p, tb1, stats] = kruskalwallis(FC);
   end
else
    fprintf('distributions are not normal\n');
   [p, tb1, stats] = kruskalwallis(FC);
end
%[p, tb1, stats] = kruskalwallis(FC);
c = multcompare(stats);
tb2 = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%only maintain connections betwwen region A and region B and between region A and
%opposite region B at the same
function [ipsi_FC, hete_FC]=maintain_specific_connections(mat)
L2R = mat(1:31,32:62);
R2L = mat(32:62,1:31);
L2L = mat(1:31,1:31);
R2R = mat(32:62,32:62);

upper_l2r = triu(L2R,1);
lower_l2r = tril(L2R, -1);
idx_l2r_u = find(upper_l2r);
idx_l2r_l = find(lower_l2r);
idx_l2r = [idx_l2r_l; idx_l2r_u];
l2l_vec = L2L(idx_l2r);
l2r_vec = L2R(idx_l2r);
if ~isempty(find(l2l_vec==0))
    fprintf('ipsi number different hete number');
    idx = find(l2l_vec);
    l2l_vec = l2l_vec(idx);
    l2r_vec = l2r_vec(idx);
end

lower_r2l = tril(R2L, -1);
upper_r2l = triu(R2L, 1);
idx_r2l_l = find(lower_r2l);
idx_r2l_u = find(upper_r2l);
idx_r2l = [idx_r2l_l; idx_r2l_u];
r2r_vec = R2R(idx_r2l);
r2l_vec = R2L(idx_r2l);
ipsi_FC =  [l2l_vec; r2r_vec];
hete_FC = [l2r_vec; r2l_vec];
if ~isempty(find(r2r_vec==0))
    fprintf('ipsi number different hete number');
    idx = find(r2r_vec);
    r2r_vec = r2r_vec(idx);
    r2l_vec = r2l_vec(idx);
end

end