function heterogeneity = HeterogeneityMutualHeterotopicDownstream(N,Nregion,W)

    heterogeneity = zeros(1,Nregion-1);
    neuron2region = zeros(Nregion,N*Nregion);
    for i = 1:Nregion
        i1 = (i-1)*N+1;
        i2 = i*N;
        neuron2region(i,:) = logical(sum(W(i1:i2,:),1));
    end
    ndownstream = zeros(1,N*Nregion);
    for i = 1:Nregion*N
        inds = 1:Nregion;
        regionid = floor(i/N)+1;
        inds = setdiff(inds,regionid);
        ndownstream(i) = sum(neuron2region(inds,i));
    end
    ndownstream = ndownstream(logical(ndownstream~=0));
    for i = 1:Nregion-1
        heterogeneity(i) = sum(logical(ndownstream==i))/length(ndownstream);
    end
end