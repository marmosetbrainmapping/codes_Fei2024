% parameters for the network.
f = 0.8;
N = 100;
Nregion = 2;

miuE = 0.1;
miuI = -miuE*f/(1-f);
sE = 0.02;
sI = sE*(f/(1-f))^0.5;
p = 0.1;
pself = 0.5;
pcontras = 0.1*ones(1,Nregion);
q = 0.1;


overlaphomo_range = 0:0.1:1;
overlaphetero = 0.5;

nE = f*N;
nI = N-nE;
ntargets = q*N;
nsources = p*N;
nself = pself*N;
ncontras = pcontras*N;
%%
Vrest = -0.1;
r = -1;
Vthres = 0;
beta = 10.;

T = 20000;
tau = 1;
dt = 0.1;
noises = [0.05,0.1,0.2,0.5];

winlen = 100;


%%
niter = 100*5;
corr_thres = 0.95;

% Cmeans: mean homotopic FC under different levels of noises and heterogeneity.
Cmeans = zeros(niter,length(noises),length(overlaphomo_range));



tic
parfor iter = 1:niter
    
    
%     tic
    disp(['### iter-',num2str(iter),', generating W...']);
    Ws = zeros(N*Nregion*2,N*Nregion*2,length(overlaphomo_range));
    for i = 1:length(overlaphomo_range)
        overlap_homos = ones(1,Nregion)*overlaphomo_range(i);
        overlap_heteros = ones(Nregion,Nregion)*overlaphetero;
        overlap_heteros(logical(eye(Nregion))) = nan;
        W = RandomNetwork_Bilateral(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,ncontras,overlap_heteros,overlap_homos);
        Ws(:,:,i) = W;
    end

%     tic
    disp(['### iter-',num2str(iter),', simulating V...']);
    Cmean = zeros(length(noises),length(overlaphomo_range));
    Cstd = Cmean;
    psync = Cmean;
    Vs = zeros(N*Nregion*2,T,length(noises),length(overlaphomo_range));
    for i = 1:length(noises)
        noise = noises(i);
        perturbation = zeros(N*Nregion*2,T);
        filtered_noise_1 = zeros(N*Nregion*2,T);
        filtered_noise_2 = zeros(N*Nregion*2,T);
        filtered_noise_3 = zeros(N*Nregion*2,T);
        for t = 1:T-1
            perturbation(:,t+1) = perturbation(:,t)+dt*(-perturbation(:,t)+filtered_noise_1(:,t));
            filtered_noise_1(:,t+1)= filtered_noise_1(:,t)+dt*0.1*(-1*filtered_noise_1(:,t)+filtered_noise_2(:,t));
            filtered_noise_2(:,t+1)= filtered_noise_2(:,t)+dt*0.1*(-1*filtered_noise_2(:,t)+filtered_noise_3(:,t));
            overall_noise = (2*noise*dt).^0.5*randn(1);
            neuron_noise = (2*noise*dt/5).^0.5*randn(N*Nregion*2,1);
            filtered_noise_3(:,t+1)= squeeze(filtered_noise_3(:,t))+squeeze(dt*(-1*filtered_noise_3(:,t)))+overall_noise+neuron_noise;
        end
        for j = 1:length(overlaphomo_range)
            W = squeeze(Ws(:,:,j));
            V = RandomNetworkDynamics(N,Nregion*2,W,Vrest,Vthres,r,beta,tau,dt,T,perturbation);
            Vregions = zeros(Nregion*2,T);
            for k = 1:Nregion*2
                r1 = (k-1)*N+1;
                r2 = k*N;
                Vregions(k,:) = mean(V(r1:r2,:),1);
            end
            C = zeros(1,T-winlen);
            for k = 1:T-winlen
                t1 = k;
                t2 = k+winlen;
                x1 = Vregions(1:Nregion,t1:t2)';
                x2 = Vregions(Nregion+1:end,t1:t2)';
                Ct = corr(x1,x2);
                Ct = Ct(logical(eye(Nregion)));
                C(k) = mean(Ct,'all','omitnan');
            end
            Cmean(i,j) = mean(abs(C(C<=corr_thres)));
            Cstd(i,j) = std(abs(C(C<=corr_thres)));
            psync(i,j) = length(find(C>=corr_thres))/(T-winlen);
        end
    end
%     toc

%     tic
    Cmeans(iter,:,:) = Cmean;
%     toc

end
toc
%%
save(['/zfs/scratch/qihang/HomoParams_Nregion',num2str(Nregion),'.mat']);



