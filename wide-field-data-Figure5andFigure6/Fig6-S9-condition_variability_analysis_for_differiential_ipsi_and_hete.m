clear;
clc;
close all;
ori_path = '/wide-field-data/new_analysis/distinguish_different_state_result/results/';
%SC_path = '/wide-field-data/structure_connectivity/full_SC_mask_wf.mat';
%SC_mask = load(SC_path).full_SC_mask_wf;
SC_path ='/wide-field-data/structure_connectivity/full_SC_type5_mask_wf.mat';
SC_mask =  load(SC_path).full_SC_type5_mask_wf;

fileID = fopen('/wide-field-data/new_analysis/distinguish_different_state_result/present_sessions.txt');
sessions = textscan(fileID, '%s');
fclose(fileID);
sessions_name = sessions{1,1};
states = cell(3,1);
states{1} = 'wake';
states{2} = 'nrem';
states{3} = 'rem';

% concatenate data by subjects order and then by states order
total_num = 2832;   
ipsi_num = 1416;  
hete_num = 1416; 
states_allsub_FC = zeros(3 * length(sessions_name), total_num);
for i=1:length(states)
    state = states{i};
    state = char(state);
    
    for j=1:length(sessions_name)
        name = sessions_name(j);
        name = char(name);
        data_path = fullfile(ori_path, name );
        file_path = dir([data_path,'/', state, '_full_FC.mat']);
        if isempty(file_path)
           fprintf('Wrong file path!\n');
           return;
        end 
       file_path = fullfile(file_path(1).folder, file_path(1).name);
       full_FC = load(file_path).full_FC;
       full_FC = full_FC .* SC_mask;
       [ipsi_FC_vector, hete_FC_vector] = maintain_specific_connections(full_FC);
       
       sample_vector = [ipsi_FC_vector; hete_FC_vector];
       m = (i -1 )*length(sessions_name) + j;
       states_allsub_FC(m,:) = sample_vector';
    end
end

conditions = zeros(3 * length(sessions_name), 3);
n = length(sessions_name);
conditions(1:n,1) = ones(length(sessions_name),1); % wake
conditions(n+1:2*n,2) = ones(length(sessions_name),1); % nrem
conditions(2*n+1:3*n,3) = ones(length(sessions_name),1);  % rem
X0 = zscore(states_allsub_FC);
Y0 = zscore(conditions);
% X0 = states_allsub_FC;
% Y0 = conditions;
n = length(sessions_name);
p = size(states_allsub_FC,2);  % independent variable
q = 3; % dependent variable
A = X0'*Y0;
[U,S,V] = svd(A);
r = min(n,min(p,q));
covariance_explained = zeros(r,1);
for i=1:r
    covariance_explained(i) = S(i,i) .^ 2;
end

covariance_sum = sum(covariance_explained);
covariance_explained = covariance_explained / covariance_sum;
covariance_explained  = cumsum(covariance_explained);
figure;
plot(1:3,covariance_explained,'marker','o');
ylim([0.5,1]);
xlabel('principle Latent variable');
ylabel('percentage of explained covariance');

sigular_values = [S(1,1),S(2,2),S(3,3)];
% bar(V);
P = permutation_test(X0,Y0,sigular_values(1),n);
disp(P);
bootstats1 = bootstrp(1000,@(x)( x), (1:n));
X0_wake = X0(1:n,:);
X0_nrem = X0(n+1:2*n,:);
X0_rem = X0(2*n+1:3*n,:);
Y0_wake = Y0(1:n,:);
Y0_nrem = Y0(n+1:2*n,:);
Y0_rem = Y0(2*n+1:3*n,:);
saliences = zeros(1000,total_num);
design_contrast = zeros(1000,3);
for a =1:1000
    index = bootstats1(a,:);
    X0_sample = [X0_wake(index,:);X0_nrem(index,:);X0_rem(index,:)];
    Y0_sample = [Y0_wake(index,:);Y0_nrem(index,:);Y0_rem(index,:)];
    B = X0_sample' * Y0_sample;
    [u,s,v] = svd(B);
    saliences(a,:) = u(:,1);
    design_contrast (a,:) = v(:,1);
end
[h, p1] = lillietest(design_contrast(:,1));
contrast_mean = mean(design_contrast,1);
contrast_sem = 2.1098* std(design_contrast,1) / sqrt(1000);  % 95% confidence interval (u-1.96*sigma/sqrt(n),u+1.96*SE)
figure;
X = categorical({'wake','nrem','rem'});
X = reordercats(X,{'wake','nrem','rem'});
b = bar(X,contrast_mean, 0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0.3 0.3 0.3];
b.CData(2,:) = [0.5 0.5 0.5];
b.CData(3,:) = [1 0.5 0];
ylabel('condition contrast');


hold on;
plot([X(1), X(1)], [contrast_mean(1) - contrast_sem(1), contrast_mean(1) + contrast_sem(1)], 'color',[0.2, 0.2, 0.2],'Linestyle', '-');
hold on;
plot([X(2), X(2)], [contrast_mean(2) - contrast_sem(2), contrast_mean(2) + contrast_sem(2)], 'color',[0.2, 0.2, 0.2],'Linestyle', '-');
hold on;
plot([X(3), X(3)], [contrast_mean(3) - contrast_sem(3), contrast_mean(3) + contrast_sem(3)], 'r','Linestyle', '-');
[p1, tb, stats1] = kruskalwallis(design_contrast);
c1 = multcompare(stats1);
saliences = abs(saliences);
saliences_ipsi = mean(saliences(:,1:ipsi_num),2);
saliences_hete = mean(saliences(:,ipsi_num+1:total_num),2);
saliences_FC = [saliences_ipsi,saliences_hete];
FC_mean = mean(saliences_FC,1);
FC_sem = std(saliences_FC,1)/ 1000;   %SE
figure;
X = categorical({'ipsi','hete'});
X = reordercats(X,{'ipsi','hete'});
b = bar(X,FC_mean, 0.5);
% ylim([0 0.05]);
b.FaceColor = 'flat';
b.CData(1,:) = [0.035 0.576 0.588];
b.CData(2,:) = [0.682 0.125 0.071];
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',7);

ylabel('condition variability(PLS saliences)');
hold on;
plot([X(1), X(1)], [FC_mean(1) - FC_sem(1), FC_mean(1) + FC_sem(1)], 'color',[0.2, 0.2, 0.2],'Linestyle', '-');
hold on;
plot([X(2), X(2)], [FC_mean(2) - FC_sem(2), FC_mean(2) + FC_sem(2)], 'color',[0.2, 0.2, 0.2],'Linestyle', '-');

% normal distribution test and variance 
[~,col]=size(saliences_FC);
h_test = zeros(col,1);
for i=1:col
    h=lillietest(saliences_FC(:,i));
    h_test(i)=h;
end
if length(find(h_test==0)) == col
   if vartestn(saliences_FC) > 0.001
       [p,tb1,stats] = anova1(FC);
   else
      fprintf('variances are not same\n');
      [p, tb1, stats] = kruskalwallis(saliences_FC);
   end
else
    fprintf('distributions are not normal\n');
   [p, tb1, stats] = kruskalwallis(saliences_FC);
end
% [p, tb1, stats] = kruskalwallis(saliences_FC);
c = multcompare(stats);
tb2 = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
function p =permutation_test(x,y,sigular_value,n)
sigular_values=zeros(1000,1);
for i=1:1000
x_order = random_permutation(x,3*n);
A = x_order'*y;
[u,s,v]= svd(A);
sigular_values(i) = s(1,1);
end
num=0;
for i=1:1000
    if sigular_values(i) > sigular_value
        num = num+1;
    end
end
p = num / 1000;

end
function data1_= random_permutation(data1,row)
index1 = randperm(row);
data1_ = data1(index1',:);

end

function [ipsi_FC, hete_FC]=maintain_specific_connections(mat)
L2R = mat(1:31,32:62);
R2L = mat(32:62,1:31);
L2L = mat(1:31,1:31);
R2R = mat(32:62,32:62);

upper_l2r = triu(L2R,1);
lower_l2r = tril(L2R,-1);
idx_l2r_1 = find(upper_l2r);
idx_l2r_2 = find(lower_l2r);
idx_l2r = [idx_l2r_1; idx_l2r_2];
l2l_vec = L2L(idx_l2r);
l2r_vec = L2R(idx_l2r);
if ~isempty(find(l2l_vec==0))
    fprintf('ipsi number different hete number');
    return;
end

upper_r2l = triu(R2L, 1);
lower_r2l = tril(R2L,-1);
idx_r2l_1 = find(upper_r2l);
idx_r2l_2 = find(lower_r2l);
idx_r2l = [idx_r2l_1; idx_r2l_2];
r2r_vec = R2R(idx_r2l);
r2l_vec = R2L(idx_r2l);
ipsi_FC =  [l2l_vec; r2r_vec];
hete_FC = [l2r_vec; r2l_vec];
if ~isempty(find(r2r_vec==0))
    fprintf('ipsi number different hete number');
    return;
end

end




