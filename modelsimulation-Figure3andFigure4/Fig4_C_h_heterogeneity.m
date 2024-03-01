f = 0.8;
N = 500;
Nregion = 30;

miuE = 0.1;
miuI = -miuE*f/(1-f);
sE = 0.02;
sI = sE*(f/(1-f))^0.5;
p = 1/Nregion;
pself = 0.5;
q = 0.1;

nE = f*N;
nI = N-nE;
ntargets = q*N;
nsources = p*N;
nself = pself*N;

Vrest = -0.1;

Vthres = 0;
beta = 10.;

T = 20000;
tau = 1;
dt = 0.1;
%%
figure
hold on
for r = 0.4:0.02:0.5

    tic
    disp(['r=',num2str(r)]);
    W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,r);
    h = HeterogeneityMutualHeterotopicDownstream(N,Nregion,W);
    toc
    plot(h,'LineWidth',2);
    

end
set(gcf,'Position',[0,0,600,400]);
set(gca,'Color','none');
set(gca,'FontSize',10);
%%
figure
hold on
for r = 0.1:0.02:0.2

    tic
    disp(['r=',num2str(r)]);
    W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,r);
    h = HeterogeneityMutualHeterotopicDownstream(N,Nregion,W);
    toc
    plot(h,'LineWidth',2);
    

end
set(gcf,'Position',[0,0,600,400]);
set(gca,'Color','none');
set(gca,'FontSize',10);
%%
figure
hold on
r1 = 0.6;
W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,r1);
h = HeterogeneityMutualHeterotopicDownstream(N,Nregion,W);
plot(h,'b','LineWidth',2);

r2 = 0.4;
W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,r2);
h = HeterogeneityMutualHeterotopicDownstream(N,Nregion,W);
plot(h,'r','LineWidth',2);

ylim([0,0.4]);
%%
figure
hold on
r1 = 0.4;
W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,r1);
h = HeterogeneityMutualHeterotopicDownstream(N,Nregion,W);
plot(h,'b','LineWidth',2);

r2 = 0.25;
W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,r2);
h = HeterogeneityMutualHeterotopicDownstream(N,Nregion,W);
plot(h,'r','LineWidth',2);


%%
r_range = [0,0.01,0.05,0.1,0.2,0.4,0.6,0.8];
hs = zeros(length(r_range),Nregion-1);
for i = 1:length(r_range)
    tic
    r = r_range(i);
    W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,r);
    h = HeterogeneityMutualHeterotopicDownstream(N,Nregion,W);
    hs(i,:) = h;
    toc
end
%%
colormap(parula(7));
cmap = colormap;
cmap = flipud(cmap);
%%
figure
hold on
legends = [];
for i = 2:length(r_range)
    h = hs(i,:);
    plot(h,'color',cmap(i-1,:),'LineWidth',1.5);
    legends = [legends,string(num2str(r_range(i)))];
end
l = legend(legends);
set(gca,'Color','none');
set(l,'Color','none');
set(gca,'FontSize',15);














