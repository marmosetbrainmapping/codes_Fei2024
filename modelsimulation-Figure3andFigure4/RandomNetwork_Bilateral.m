function [W,projectingneurons] = RandomNetwork_Bilateral(N,Nregion,miuE,sE,nE,miuI,sI,ntargets,nsources,nself,ncontras,overlap_heteros,overlap_homos)

%   N: the number of neurons per region.
%   Nregion: the number of regions.
%   miuE(sE): mean(s.t.d) excitatory synaptic strength.
%   miuI(sI): mean(s.t.d) inhibitory synaptic strength.
%   ntargets: the number of target neurons per region for each upstream
% projecting neuron
%   nsources: the number of upstream projecting neurons for each pair of
% upstream-downstream region.
%   nself: the number of neurons per region projecting to the region itself.
%   ncontras (1xNregion): the number of neurons projecting to contralateral homotopic
% region.
%   overlap_heteros (NregionxNregion): overlap of bilateral heterotopic projections.
%   overlap_homos (1xNregion): overlap of self-projection and homotopic projection.

    W = zeros(N*Nregion*2,N*Nregion*2);
    projectingneurons = cell(Nregion*2,Nregion*2);
    
    for i = 1:Nregion*2
    
        Wi = zeros(N*Nregion*2,N);
        if i <= Nregion
            contra_i = i+Nregion;
            ncontra = ncontras(i);
            overlap_homo = overlap_homos(i);
        else
            contra_i = i-Nregion;
            ncontra = ncontras(contra_i);
            overlap_homo = overlap_homos(contra_i);
        end

        randinds = randperm(N);
        Einds = randinds(1:nE);
        
        sourceinds = cell(Nregion*2,1);
        for j = 1:Nregion*2
            if j == i
                % generate self-projecting neurons (selfinds)
                selfcontra_inds = randperm(N);
                selfinds = selfcontra_inds(1:nself);
                sourceinds{j} = selfinds;
                % generate homotopic-projecting neurons (contrainds)
                if overlap_homo ~= -1
                    ind1 = nself-round(ncontra*overlap_homo)+1;
                    ind2 = ind1+ncontra-1;
                    contrainds = selfcontra_inds(ind1:ind2);
                    sourceinds{contra_i} = contrainds;
                else
                    contrainds = randperm(N);
                    contrainds = contrainds(1:ncontra);
                    sourceinds{contra_i} = contrainds;
                end
            elseif j == contra_i
                continue
            elseif j <= Nregion
                % generate heterotopic-projecting neurons
                if i <= Nregion
                    overlap_hetero = overlap_heteros(j,i);
                else
                    overlap_hetero = overlap_heteros(j,contra_i);
                end
                contra_j = j+Nregion;
                heterotopic_inds = randperm(N);
                sourceinds{j} = heterotopic_inds(1:nsources);
                ind1 = nsources-round(nsources*overlap_hetero)+1;
                ind2 = ind1+nsources-1;
                sourceinds{contra_j} = heterotopic_inds(ind1:ind2);
            else
                continue
            end
        end
        for j = 1:Nregion*2
            projectingneurons{j,i} = sourceinds{j};
        end

    
        for j = 1:Nregion*2
            Wij = zeros(N,N);
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
           
    
            Wi((j-1)*N+1:j*N,:) = Wij;
        end
    
        W(:,(i-1)*N+1:i*N) = Wi;
    
    end
    
    for i = 1:N*Nregion*2
        input = W(i,:);
        inds = find(input~=0);
        W(i,inds) = W(i,inds) - sum(W(i,:))/length(inds);
    end

end