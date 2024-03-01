clear all
%%
f = 0.8;
N = 100;
Nregion_range = [2,3,5,7,10];



miuE = 0.1;
miuI = -miuE*f/(1-f);
sE = 0.02;
sI = sE*(f/(1-f))^0.5;
p = 0.05;
pself = 0.5;
q = 0.05;


overlaphomo = 0.5;
overlaphetero = 0.5;
pcontra1 = 0;
pcontra2 = 1;

nE = f*N;
nI = N-nE;
ntargets = q*N;
nsources = p*N;
nself = pself*N;
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
niter = 100*3;
corr_thres = 0.95;

variation_range = zeros(niter,length(noises),length(Nregion_range));


tic
parfor iter = 1:niter
    
    
%     tic
    disp(['### iter-',num2str(iter),', generating W...']);

    Ws = cell(length(Nregion_range),2);
    for i = 1:length(Nregion_range)
        Nregion = Nregion_range(i);
        pcontras = pcontra1*ones(1,Nregion);
        ncontras = pcontras*N;
        overlap_homos = ones(1,Nregion)*overlaphomo;
        overlap_heteros = ones(Nregion,Nregion)*overlaphetero;
        overlap_heteros(logical(eye(Nregion))) = nan;
        W = RandomNetwork_Bilateral(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,ncontras,overlap_heteros,overlap_homos);
        Ws{i,1} = W;

        pcontras = pcontra2*ones(1,Nregion);
        ncontras = pcontras*N;
        W = RandomNetwork_Bilateral(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,ncontras,overlap_heteros,overlap_homos);
        Ws{i,2} = W;
    end

%     tic
    disp(['### iter-',num2str(iter),', simulating V...']);
    vr = zeros(length(noises),length(Nregion_range));
    for i = 1:length(noises)
        noise = noises(i);
        
        for j = 1:length(Nregion_range)
            Nregion = Nregion_range(j);
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

            W = Ws{j,1};
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
            C1 = mean(abs(C(C<=corr_thres)));

            W = Ws{j,2};
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
            C2 = mean(abs(C(C<=corr_thres)));

            vr(i,j) = C2-C1;
        end
    end
%     toc

%     tic
    variation_range(iter,:,:) = vr;
%     toc

end
toc
%%
save('/zfs/scratch/qihang/PcontrasVariationRange.mat');


