function W = RandomNetwork(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,ratio)

    W = zeros(N*Nregion,N*Nregion);
    nsources_min = nsources;
    nsources_max = (Nregion-1)*nsources;
    nsources_total = nsources_min + (nsources_max-nsources_min)*ratio;
    
    for i = 1:Nregion
    
        Wi = zeros(N*Nregion,N);
        randinds = randperm(N);
        Einds = randinds(1:nE);
        
        randinds = randperm(N);
        hetero_sources = randinds(1:nsources_total);
        sourceinds = cell(Nregion,1);
        for j = 1:Nregion
            sourceinds{j} = [];
        end
        selfinds = randperm(N);
        selfinds = selfinds(1:nself);
        sourceinds{i} = selfinds;
        
        hetero_regions = setdiff((1:Nregion),i);
        for j = 1:nsources_total
            neuron = hetero_sources(j);
            regionidx = randperm(length(hetero_regions));
            region = hetero_regions(regionidx(1));
            sourceinds{region} = [sourceinds{region},neuron];
            fullregions = [];
            for k = 1:length(hetero_regions)
                if length(sourceinds{hetero_regions(k)}) == nsources
                    fullregions = [fullregions,hetero_regions(k)];
                end
            end
            hetero_regions = setdiff(hetero_regions,fullregions);
        end
        
        for j = 1:Nregion
            if j == i
                continue
            end
            heteroinds = sourceinds{j};
            tmpsources = setdiff(hetero_sources,heteroinds);
            nlack = nsources - length(heteroinds);
            randinds = randperm(length(tmpsources));
            randinds = randinds(1:nlack);
            heteroinds = [heteroinds,tmpsources(randinds)];
            sourceinds{j} = heteroinds;
        end
        
%         for j = 1:Nregion
%             if j == i
%                 selfinds = randperm(N);
%                 selfinds = selfinds(1:nself);
% %                 sourceinds(j,:) = selfinds;
%                 sourceinds{j} = selfinds;
%             else
%                 if ratio ~= -1
%                     ind1 = round(count*nsources)+1;
%                     ind2 = round(count*nsources)+nsources;
%                     otherinds = randinds(ind1:ind2);
%                     count = count + ratio;
%                     sourceinds{j} = otherinds;
%                 elseif ratio == -1
%                     otherinds = randperm(N);
%                     otherinds = otherinds(1:nsources);
%                     sourceinds{j} = otherinds;
%                 end
%             end
%         end
    
        for j = 1:Nregion
            Wij = zeros(N,N);
%             inds = sourceinds(j,:);
            inds = sourceinds{j};
            synapses = zeros(N,length(inds));
            parfor k = 1:length(inds)
                ind = inds(k);
                targetinds = randperm(N);
                targetinds = targetinds(1:ntargets);
                if ismember(ind,Einds)
                    synapse = zeros(N,1);
                    synapse(targetinds) = randn(ntargets,1)*sE + miuE;
                    synapses(:,k) = synapse;
                else
                    synapse = zeros(N,1);
                    synapse(targetinds) = randn(ntargets,1)*sI + miuI;
                    synapses(:,k) = synapse;
                end
            end
            Wij(:,inds) = synapses;
            if j == i
                Wij(logical(eye(N,N))) = 0;
            end
            
            

%             for k = inds
%                 targetinds = randperm(N);
%                 targetinds = targetinds(1:ntargets);
%                 if ismember(k,Einds)
%                     synapse = randn(ntargets,1)*sE + miuE;
%                     Wij(targetinds,k) = synapse;
%                     Wij(logical(eye(N,N))) = 0;
%                 else
%                     synapse = randn(ntargets,1)*sI + miuI;
%                     Wij(targetinds,k) = synapse;
%                     Wij(logical(eye(N,N))) = 0;
%                 end
%             end
    
            Wi((j-1)*N+1:j*N,:) = Wij;
        end
    
        W(:,(i-1)*N+1:i*N) = Wi;
    
    end
    
    for i = 1:N*Nregion
        input = W(i,:);
        inds = find(input~=0);
        W(i,inds) = W(i,inds) - sum(W(i,:))/length(inds);
    end

end