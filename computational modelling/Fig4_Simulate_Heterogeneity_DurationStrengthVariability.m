% parameters for the network.
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
noises = [0.01,0.02,0.05,0.1,0.2,0.5];
ratios = 0:0.1:1;
%%
niter = 1000;
winlen = 100;
corr_thres = 0.90;


% Cmeans: mean inter-regional FC over time under different levels of noises and heterogeneity.
% Cstds: inter-regional FC standard deviation over time.
% psyncs: duration of synchronous periods.
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
