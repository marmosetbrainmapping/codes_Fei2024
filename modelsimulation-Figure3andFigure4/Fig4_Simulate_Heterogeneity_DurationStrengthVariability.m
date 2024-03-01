f = 0.8;
N = 100;
Nregion = 5;

miuE = 0.1;
miuI = -miuE*f/(1-f);
sE = 0.02;
sI = sE*(f/(1-f))^0.5;
p = 0.25;
pself = 0.5;
q = 0.2;

nE = f*N;
nI = N-nE;
ntargets = q*N;
nsources = p*N;
nself = pself*N;

Vrest = -0.1;
r = -1;
Vthres = 0;
beta = 10.;

T = 20000;
tau = 1;
dt = 0.1;
%%
% ratios = [0,0.25,0.5,0.75,1];
% W_ratios = zeros(N*Nregion,N*Nregion,length(ratios));
% 
% for i = 1:length(ratios)
%     ratio = ratios(i);
%     W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,ratio);
%     W_ratios(:,:,i) = W;
% end
%%
noises = [0.01,0.02,0.05,0.1,0.2,0.5];
ratios = 0:0.1:1;
%%
niter = 1000;
winlen = 100;
corr_thres = 0.90;

Cmeans = zeros(niter,length(noises),length(ratios));
Cstds = zeros(niter,length(noises),length(ratios));
psyncs = zeros(niter,length(noises),length(ratios));

tic
parfor iter = 1:niter
    
    
%     tic
    disp(['### iter-',num2str(iter),', generating W...']);
    Ws = zeros(N*Nregion,N*Nregion,length(ratios));
    for i = 1:length(ratios)
        ratio = ratios(i);
        W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,ratio);
        Ws(:,:,i) = W;
    end
%     toc

    tic
    disp(['### iter-',num2str(iter),', simulating V...']);
    Vs = zeros(N*Nregion,T,length(noises),length(ratios));
    for i = 1:length(noises)
        noise = noises(i);
        perturbation = zeros(N*Nregion,T);
        filtered_noise_1 = zeros(N*Nregion,T);
        filtered_noise_2 = zeros(N*Nregion,T);
        filtered_noise_3 = zeros(N*Nregion,T);
        for t = 1:T-1
            perturbation(:,t+1) = perturbation(:,t)+dt*(-perturbation(:,t)+filtered_noise_1(:,t));
            filtered_noise_1(:,t+1)= filtered_noise_1(:,t)+dt*0.1*(-1*filtered_noise_1(:,t)+filtered_noise_2(:,t));
            filtered_noise_2(:,t+1)= filtered_noise_2(:,t)+dt*0.1*(-1*filtered_noise_2(:,t)+filtered_noise_3(:,t));
            overall_noise = (2*noise*dt).^0.5*randn(1);
            neuron_noise = (2*noise*dt/5).^0.5*randn(N*Nregion,1);
            filtered_noise_3(:,t+1)= squeeze(filtered_noise_3(:,t))+squeeze(dt*(-1*filtered_noise_3(:,t)))+overall_noise+neuron_noise;
        end
        for j = 1:length(ratios)
            W = squeeze(Ws(:,:,j));
            V = RandomNetworkDynamics(N,Nregion,W,Vrest,Vthres,r,beta,tau,dt,T,perturbation);
            Vs(:,:,i,j) = V;
        end
    end
    toc

%     tic
    disp(['### iter-',num2str(iter),', calculating C...']);
    Cmean = zeros(length(noises),length(ratios));
    Cstd = Cmean;
    psync = Cmean;
    for i = 1:length(noises)
        for j = 1:length(ratios)
            V = squeeze(Vs(:,:,i,j));
            Vregions = zeros(Nregion,T);
            for k = 1:Nregion
                r1 = (k-1)*N+1;
                r2 = k*N;
                Vregions(k,:) = mean(V(r1:r2,:),1);
            end
            C = zeros(1,T-winlen);
            for k = 1:T-winlen
                t1 = k;
                t2 = k+winlen;
                x = Vregions(1:Nregion,t1:t2)';
                Ct = corr(x);
                Ct = Ct(~logical(eye(Nregion)));
                C(k) = mean(Ct,'all','omitnan');
            end
            Cmean(i,j) = mean(abs(C(C<=corr_thres)));
            Cstd(i,j) = std(abs(C(C<=corr_thres)));
            psync(i,j) = length(find(C>=corr_thres))/(T-winlen);
        end
    end
    Cmeans(iter,:,:) = Cmean;
    Cstds(iter,:,:) = Cstd;
    psyncs(iter,:,:) = psync;
%     toc

end
toc

save(['/zfs/scratch/qihang/ratioParams_thres',num2str(corr_thres),'.mat']);
%%
load(['/zfs/scratch/qihang/ratioParams_thres',num2str(corr_thres),'.mat']);
%%
markers = ["o","*","diamond","^"];
lines = ["-","--",":","-."];
sz = 50;
left_color = [0.4940 0.1840 0.5560];
right_color = [0.4660 0.6740 0.1880];
fig = figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on

for i = 1:length(noises)

yyaxis left
x = ratios;
y = mean(abs(Cmeans(:,i,:)),1);
err = std(abs(Cmeans(:,i,:)),[],1)/size(Cmeans(:,i,:),1)^0.5;
errorbar(x,squeeze(y),squeeze(err),'Color',left_color,'LineWidth',1.5,'LineStyle',lines(i));
s = scatter(x,squeeze(y),sz,markers(i));
% s.MarkerFaceColor = none;
s.MarkerEdgeColor = left_color;

yyaxis right
x = ratios;
y = mean(abs(Cstds(:,i,:)),1);
err = std(abs(Cstds(:,i,:)),[],1)/size(Cmeans(:,i,:),1)^0.5;
errorbar(x,squeeze(y),squeeze(err),'Color',right_color,'LineWidth',1.5,'LineStyle',lines(i));
s = scatter(x,squeeze(y),sz,markers(i));
% s.MarkerFaceColor = none;
s.MarkerEdgeColor = right_color;

l = legend(["mean","s.t.d"]);
xlabel('divergence');
axis square
axis tight
set(gca,'Color','none')
set(l,'Color','none')
% set(gcf,'Color','none')
end
%%


load(['/zfs/scratch/qihang/ratioParams_thres',num2str(corr_thres),'.mat']);
% colormap(parula(length(noises)-1));
% cmap = colormap;
% cmap = flipud(cmap);


%%
lines = ["-","--",":","-.","--",":"];
markers = ["o","*","diamond","^",">","<"];
left_colors = ["#7E2F8E","#904D9E","#A26BAE","#B489BE","#C6A7CE","#D8C5DE"];
right_colors = ["#77AC30","#8AB84E","#9DC46C","#B0D08A","#C3DCA8","#D6E8C6"];
sz = 8;
fig = figure;
set(fig,'defaultAxesColorOrder',[[126,47,142]/255; [119,172,48]/255]);
hold on
for i = 1:length(noises)-1
    yyaxis left
    x = ratios;
    y = mean(abs(Cmeans(:,i,:)),1);
    plot(x,squeeze(y)/max(y),'Color',left_colors(i),'LineWidth',1.5,'MarkerEdgeColor',left_colors(i),'Marker',markers(i),'MarkerSize',sz,'LineStyle',lines(i));


    yyaxis right
    x = ratios;
    y = mean(abs(Cstds(:,i,:)),1);
    plot(x,squeeze(y)/max(y),'Color',right_colors(i),'LineWidth',1.5,'MarkerEdgeColor',right_colors(i),'Marker',markers(i),'MarkerSize',sz,'LineStyle',lines(i));

end
%%
figure
hold on
x = ratios;
ys = squeeze(mean(abs(psyncs),1));
legends = cell(length(noises),1);
for i = 1:length(noises)
    y = ys(i,:);
    plot(x,y/max(y),'LineWidth',2);
    legends{i} = num2str(noises(i));
end
legend(legends)
%%
lines = ["-","--",":","-.","--",":"];
markers = ["o","*","diamond","^",">","<"];
left_colors = ["#7E2F8E","#904D9E","#A26BAE","#B489BE","#C6A7CE","#D8C5DE"];
right_colors = ["#77AC30","#8AB84E","#9DC46C","#B0D08A","#C3DCA8","#D6E8C6"];
sz = 8;
fig = figure;
set(fig,'defaultAxesColorOrder',[[126,47,142]/255; [119,172,48]/255]);
hold on
i = 3;
yyaxis left
x = ratios;
y = mean(abs(Cmeans(:,i,:)),1);
plot(x,squeeze(y)/max(y),'Color',left_colors(1),'LineWidth',1.5,'MarkerEdgeColor',left_colors(1),'Marker',markers(1),'MarkerSize',sz,'LineStyle',lines(1));


yyaxis right
x = ratios;
y = mean(abs(Cstds(:,i,:)),1);
plot(x,squeeze(y)/max(y),'Color',right_colors(1),'LineWidth',1.5,'MarkerEdgeColor',right_colors(1),'Marker',markers(1),'MarkerSize',sz,'LineStyle',lines(1));
%%
load('/zfs/scratch/qihang/ratioParams.mat');
% colormap(turbo(length(noises)-1+2));
% cmap = colormap;
% cmap = flipud(cmap);

cmap = colormap(cbrewer2('RdBu', length(noises)-1));
cmap = flipud(cmap);

lines = ["-","--",":","-.","--",":"];
markers = ["o","*","diamond","^",">","<"];
sz = 8;
%%
figure
hold on
for i = 1:length(noises)-1

    x = ratios;
    y = mean(abs(Cmeans(:,i,:)),1);
    plot(x,squeeze(y)/max(y),'Color',cmap(i,:),'LineWidth',1.5,'MarkerEdgeColor',cmap(i,:),'Marker',markers(i),'MarkerSize',sz,'LineStyle',lines(i),'MarkerFaceColor','none');

end
set(gca,'Color','none');
set(gca,'FontSize',15);
set(gcf,'Position',[0,0,600,400]);
set(gcf,'Color','none');
%%
figure
hold on
for i = 1:length(noises)-1

    x = ratios;
    y = mean(abs(Cstds(:,i,:)),1);
    plot(x,squeeze(y)/max(y),'Color',cmap(i,:),'LineWidth',1.5,'MarkerEdgeColor',cmap(i,:),'Marker',markers(i),'MarkerSize',sz,'LineStyle',lines(i),'MarkerFaceColor','none');

end
set(gca,'Color','none');
set(gca,'FontSize',15);
set(gcf,'Position',[0,0,600,400]);
set(gcf,'Color','none');
%%
figure
hold on
legends = [];
for i = 1:length(noises)-1

    x = ratios;
    y = mean(abs(psyncs(:,i,:)),1);
    plot(x,squeeze(y)/max(y),'Color',cmap(i,:),'LineWidth',1.5,'MarkerEdgeColor',cmap(i,:),'Marker',markers(i),'MarkerSize',sz,'LineStyle',lines(i),'MarkerFaceColor','none');
    legends = [legends,string(num2str(noises(i)))];

end
l = legend(legends);
set(l,'Color','none');
set(gca,'Color','none');
set(gca,'FontSize',15);
set(gcf,'Position',[0,0,600,400]);
set(gcf,'Color','none');
%%













