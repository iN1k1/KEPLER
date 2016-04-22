function [ cmc, finalDist, cmcFeat, featDist, ...
    score, scoreFeat ] = my_calcMCMC( model, data, ...
    testSet, weights, fusionFunction, onlyScores, multipleShotFusionFunction )

if nargin < 5
    fusionFunction = @sum;
    onlyScores = false;
elseif nargin < 6
    onlyScores = false;
end
    
idxatest = unique(testSet.index(:,1), 'stable');
idxbtest = unique(testSet.index(:,2), 'stable');

featDist = zeros([length(idxatest) length(idxbtest) length(model)]);
weightedDistFeat = featDist;
distToRemove = [];
parfor f=1:length(model)
    
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
    
    M = model(f).M;
    featDist(:,:,f) = logistic(0.02 * sqdist(Xa, Xb,M));
    weightedDistFeat(:,:,f) = weight .* featDist(:,:,f);
   
    if all(isnan(weight))
        distToRemove = [distToRemove; f];
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
    cmc = NM_CMC_ROC(IDs, score, 'ROC', false, 'isDist', true, 'statsFunHandle', multipleShotFusionFunction); 
end

% Individual features score
featD = squeeze(featDist(:,:,end));
scoreFeat = zeros(size(featD(inds),1), size(featDist,3));
cmcFeat = zeros(length(cmc), size(featDist,3));
parfor f=1:size(featDist,3)
    featD = squeeze(featDist(:,:,f));
    scoreFeat(:,f) = featD(inds);
    if ~onlyScores
        cmcFeat(:,f) = NM_CMC_ROC(IDs, scoreFeat(:,f), 'ROC', false, 'isDist', true, 'statsFunHandle', multipleShotFusionFunction); 
    end
end

end
