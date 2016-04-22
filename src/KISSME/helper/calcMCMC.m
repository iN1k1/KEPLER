function [ cmc, finalDist, cmcFeat, featDist, ...
    score, scoreFeat ] = calcMCMC( model, data, ...
    testSet, weights, fusionFunction, onlyScores, multipleShotFusionFunction )

if nargin < 5
    fusionFunction = @sum;
    onlyScores = false;
elseif nargin < 6
    onlyScores = false;
end
    
%idatest = unique(testSet.ID(:,1), 'stable');
%idbtest = unique(testSet.ID(:,2), 'stable');
idxatest = unique(testSet.index(:,1), 'stable');
idxbtest = unique(testSet.index(:,2), 'stable');
%idxapos = arrayfun(@(idxa)(find(idxa==idxatest)), testSet.index(:,1));
%idxbpos = arrayfun(@(idxb)(find(idxb==idxbtest)), testSet.index(:,2));

featDist = zeros([length(idxatest) length(idxbtest) length(model)]);
weightedDistFeat = featDist;
distToRemove = [];
for f=1:length(model)
    
    if iscell(data)
        if isempty(data{f})
            continue;
        end
    
        Xa = data{f}(1:size(model(f).M,1),idxatest);
        Xb = data{f}(1:size(model(f).M,1),idxbtest);
    else
        Xa = data(1:size(model(f).M,1),idxatest);
        Xb = data(1:size(model(f).M,1),idxbtest);
    end
    
    if nargin >= 3 && ~isempty(weights)
        weight = weights(f);
    else
        weight = 1/length(model);
    end
    
    % 
    % if isfield(s, 'A')
    %     Xa = Xa';
    %     Xb = Xb';
    %     N = size(Xa,1);
    %     Xa = (Xa-repmat(s.muA,N,1))*s.A;
    %     Xb = (Xb-repmat(s.muB,N,1))*s.B;
    %     Xa = Xa';
    %     Xb = Xb';
    % end
    if isfield(model, 'mapping')
        Xa = out_of_sample(Xa', model.mapping);
        Xb = out_of_sample(Xb', model.mapping);
        Xa = Xa';
        Xb = Xb';
    end

%     if isfield(s, 'Wa')
%         Xa = bsxfun(@minus, Xa, s.muA);
%         Xb = bsxfun(@minus, Xb, s.muB);
% 
%         Xa = Xa'*s.Wa;
%         Xb = Xb'*s.Wb;
% 
%         Xa = Xa';
%         Xb = Xb';
% 
%         distFeat = NM_pdist(Xa', Xb', 'chisq');
%     end

    if isfield(model, 'gmm')

        Xdiff = zeros(size(Xa,2)*size(Xb,2), size(Xa,1));
        pairDist = zeros(size(Xa,2)*size(Xb,2), 1);
        pairIDs = zeros(size(Xa,2)*size(Xb,2), 2);
        pariABidx = zeros(size(Xa,2)*size(Xb,2), 2);
        i = 1;
        for a=1: size(Xa,2)
            for b=1:size(Xb,2)
                xa = Xa(:,a);
                xb = Xb(:,b);
                Xdiff(i,:) = abs(xa-xb);
                pairIDs(i,:) = [idxatest(a) idxbtest(b)];
                pariABidx(i,:) = [a b];
                i = i + 1;
            end
        end

        featDist = zeros(size(Xa,2), size(Xb,2));
        [idx,nlogl] = cluster(model.gmm,Xdiff);

        for ab=1:length(pariABidx)

                xa = Xa(:,pariABidx(ab,1));
                xb = Xb(:,pariABidx(ab,2));

                posteriorProb = model.gmm.posterior(Xdiff(ab,:));

                for i=1:model.gmm.NComponents
        %         pos = find(idx==i);
        %         idxInCurrentCluster = pairIDs(pos,:);
        %         abIdxInCurrentCluster = pariABidx(pos,:);
        %         for ab=1:length(idxInCurrentCluster)
        %             xa = Xa(:,abIdxInCurrentCluster(ab,1));
        %             xb = Xb(:,abIdxInCurrentCluster(ab,2));
        %             
        %             pairDist(pos(ab)) = NM_pdist(xa',xb', 'quaddiff', s.M{i});
        %             
        %             dist(abIdxInCurrentCluster(ab,1), abIdxInCurrentCluster(ab,2)) = pairDist(pos(ab));
        %         end

                    featDist(pariABidx(ab,1), pariABidx(ab,2)) = ...
                        featDist(pariABidx(ab,1), pariABidx(ab,2))...
                        + (NM_pdist(xa',xb', 'quaddiff', model.M{i}) * posteriorProb(i));
                end

        end

    %    pars.dataset.name = 'VIPeR';
    %    result = NM_CMC_ROC( pairIDs, ones(size(Xdiff,1),1), max(pairDist)-pairDist, pars, 'ROC', false);
    elseif isfield(model, 'MT')
        Xte = model(f).normData;
        idxtestall = [];
        for a=1:length(idxatest)
            for b=1:length(idxbtest)
                scoreTT(a,b) = Score_SML(Xte, idxatest(a), idxbtest(b), model(f).MT, model(f).GT);
            end
        end
        featDist(:,:,f) = max(scoreTT(:)) - scoreTT;     
        featDist(:,:,f) = weight * featDist(:,:,f);
        
        %ScoreTTPOS = Score_SML(Xte, SS_te(:,1), SS_te(:,2), s.MT, s.GT);
        %ScoreTTNEG = Score_SML(Xte, DD_te(:,1), DD_te(:,2), s.MT, s.GT);
        %ScoreTT = [ScoreTTPOS; ScoreTTNEG];
    else
        M = model(f).M;
        featDist(:,:,f) = logistic(0.02 * sqdist(Xa, Xb,M));
        weightedDistFeat(:,:,f) = weight .* featDist(:,:,f);
        %dist = dist + NM_pdist(Xa',Xb', 'cosine');
        
        if all(isnan(weight))
            distToRemove = [distToRemove; f];
        end
    end
end

% Remove distance not used for re-id
weightedDistFeat(:,:,distToRemove) = [];

% Compute final distance
finalDist = fusionFunction(weightedDistFeat,3);
idxapos = arrayfun(@(idxa)(find(idxa==idxatest)), testSet.index(:,1));
idxbpos = arrayfun(@(idxb)(find(idxb==idxbtest)), testSet.index(:,2));
inds = sub2ind(size(finalDist), idxapos, idxbpos);
score = finalDist(inds);

% Compute CMC
IDs = testSet.ID;
cmc = zeros(length(unique(IDs(:,2))),1);
if ~onlyScores
    cmc = NM_CMC_ROC(IDs, IDs, score, '', false, true, false, multipleShotFusionFunction); 
end

% Individual features score
featD = squeeze(featDist(:,:,f));
scoreFeat = zeros(size(featD(inds),1), size(featDist,3));
cmcFeat = zeros(length(cmc), size(featDist,3));
parfor f=1:size(featDist,3)
    featD = squeeze(featDist(:,:,f));
    scoreFeat(:,f) = featD(inds);
    if ~onlyScores
        cmcFeat(:,f) = NM_CMC_ROC(IDs, IDs, scoreFeat(:,f), '', false, true, false, multipleShotFusionFunction); 
    end
end




% testSet.scoreFeat = zeros(size(testSet(i,t).ID,1),2,size(featDists{i,t},3));
% for f=1:size(featDists{i,t},3)
%     %idxapos = arrayfun(@(idxa)(find(idxa==idxatest)), testSet(i,t).index(:,1));
%     %idxbpos = arrayfun(@(idxb)(find(idxb==idxbtest)), testSet(i,t).index(:,2));
%     featD = squeeze(featDists{i,t}(:,:,f));
%     inds = sub2ind(size(featD), idxapos, idxbpos);
%     tests(i,t).scoreFeat(:,2,f) = featD(inds);
% end
% for pairCounter=1:size(finalDist,1)
%     distPair = finalDist(pairCounter,:);  
%     [~,idx] = sort(distPair,'ascend');
%     matchRank = idatest(pairCounter) == idbtest(idx);
%     cmc(matchRank) = cmc(matchRank) + 1;
%     
%     for f=1:length(model)
%         distPair = squeeze(distFeat(pairCounter,:,f)); 
%         [tmp,idx] = sort(distPair,'ascend');
%         matchRank = idatest(pairCounter) == idbtest(idx);
%         cmcFeat(1,matchRank,f) = cmcFeat(1,matchRank,f) + 1;
%     end
%     
% end
%     
% for f=1:length(model)
%     cmcFeat(1,:,f) = cumsum(cmcFeat(1,:,f));
% end
%     
% cmc = cumsum(cmc);
% cmc = (cmc./max(cmc))*100;
% for f=1:length(model)
%     c = cmcFeat(1,:,f);
%     cmcFeat(1,:,f) = (c / max(c(:))) * 100;
% end

end
