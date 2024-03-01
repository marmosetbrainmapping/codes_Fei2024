f = 0.8;
N = 100;
Nregion = 5;

miuE = 0.1;
miuI = -miuE*f/(1-f);
sE = 0.02;
sI = sE*(f/(1-f))^0.5;
p = 0.25;
pself = 0.5;
q = 0.1;

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

Noise = 0.1;
perturbation = zeros(N*Nregion,T);
filtered_noise_1 = zeros(N*Nregion,T);
filtered_noise_2 = zeros(N*Nregion,T);
filtered_noise_3 = zeros(N*Nregion,T);
for t = 1:T-1
    perturbation(:,t+1) = perturbation(:,t)+dt*(-perturbation(:,t)+filtered_noise_1(:,t));
    filtered_noise_1(:,t+1)= filtered_noise_1(:,t)+dt*0.1*(-1*filtered_noise_1(:,t)+filtered_noise_2(:,t));
    filtered_noise_2(:,t+1)= filtered_noise_2(:,t)+dt*0.1*(-1*filtered_noise_2(:,t)+filtered_noise_3(:,t));
    overall_noise = (2*Noise*dt).^0.5*randn(1);
    neuron_noise = (2*Noise*dt*0.01).^0.5*randn(N*Nregion,1);
    filtered_noise_3(:,t+1)= squeeze(filtered_noise_3(:,t))+squeeze(dt*(-1*filtered_noise_3(:,t)))+overall_noise+neuron_noise;
end
%%
Whomo = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,0);
Wmid = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,0.5);
Whetero = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,1);
%%
Vhomo = RandomNetworkDynamics(N,Nregion,Whomo,Vrest,Vthres,r,beta,tau,dt,T,perturbation);
Vmid = RandomNetworkDynamics(N,Nregion,Wmid,Vrest,Vthres,r,beta,tau,dt,T,perturbation);
Vhetero = RandomNetworkDynamics(N,Nregion,Whetero,Vrest,Vthres,r,beta,tau,dt,T,perturbation);
%%
ymax = 2;
t1 = 5000;
t2 = 15000;

%%
figure
hold on
V = Vhomo;
for i = 1:50:N*Nregion
    plot(V(i,t1:t2),'color',[0.,0.,0.],'LineStyle','-','LineWidth',0.5);
end
plot(mean(V(:,t1:t2),1),'r','LineStyle','-','LineWidth',2);
set(gcf,'Position',[0         684        1000         150]);
set(gca,'Position',[0.02,0.1,0.95,0.8]);
title('heterogeneity = 0');
ylim([-ymax,ymax]);
xlim([1,t2-t1]);
set(gca,'Color','none');
%%
figure
hold on
V = Vmid;
for i = 1:50:N*Nregion
    plot(V(i,t1:t2),'color',[0.,0.,0.],'LineStyle','-','LineWidth',0.5);
end
plot(mean(V(:,t1:t2),1),'r','LineStyle','-','LineWidth',2);
set(gcf,'Position',[0         684        1000         150]);
set(gca,'Position',[0.02,0.1,0.95,0.8]);
title('heterogeneity = 0.5');
ylim([-ymax,ymax]);
xlim([1,t2-t1]);
set(gca,'Color','none');
%%
figure
hold on
V = Vhetero;
for i = 1:50:N*Nregion
    plot(V(i,t1:t2),'color',[0.,0.,0.],'LineStyle','-','LineWidth',0.5);
end
plot(mean(V(:,t1:t2),1),'r','LineStyle','-','LineWidth',2);
set(gcf,'Position',[0         684        1000         150]);
set(gca,'Position',[0.02,0.1,0.95,0.8]);
title('heterogeneity = 1');
ylim([-ymax,ymax]);
xlim([1,t2-t1]);
set(gca,'Color','none');
%%
W = Whomo;
h = HeterogeneityMutualHeterotopicDownstream(100,5,W);
figure
plot(h,'black','LineWidth',2);
set(gcf,'Position',[0,0,300,200]);
set(gca,'Color','none');
set(gca,'FontSize',10);
set(gcf,'Color','none');

%%

W = Wmid;
h = HeterogeneityMutualHeterotopicDownstream(100,5,W);
figure
plot(h,'black','LineWidth',2);
set(gcf,'Position',[0,0,300,200]);
set(gca,'Color','none');
set(gca,'FontSize',10);
set(gcf,'Color','none');
%%
W = Whetero;
h = HeterogeneityMutualHeterotopicDownstream(100,5,W);
figure
plot(h,'black','LineWidth',2);
set(gcf,'Position',[0,0,300,200]);
set(gca,'Color','none');
set(gca,'FontSize',10);
set(gcf,'Color','none');


%%






