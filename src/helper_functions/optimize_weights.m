function [weights] = optimize_weights(X, dataset, pars, model, trainSet)

% If all weights are nan => go for no MML
if ~isempty(pars.classifier.weights) && all(isnan(pars.classifier.weights))
    nFeats = length(model(1).kissme);
    weights = nan(size(pars.settings.testCams,1), pars.settings.numTests, nFeats);
    weights(:,:,end) = 1;
    return;
end

% Evaluate models with training data
testPars = pars;
testPars.classifier.weights = ones(size(pars.settings.testCams,1), length(X)) / length(X);
[tests, finalDists, cmcs, featDists, cmcFeats] = test_dist(model, X, trainSet, dataset, testPars, false, false);

% Exclude concatenated features
nPars = size(featDists{1},3)-1;
rng(2);

% Outputs
weights = zeros(size(pars.settings.testCams,1), pars.settings.numTests, nPars);

if nPars == 1
    % Exclude concatnated features
    weights(:) = 1;
    weights = cat(3, weights, nan(size(weights,1), size(weights,2)));
    return;
end

for i=1:size(pars.settings.testCams,1)
    for t=1:pars.settings.numTests

        % Compute Rank 1
        cmcs = cmcFeats{i,t};
        r1 = zeros(1,size(cmcs,2));
        exps = zeros(1,size(cmcs,2));
        for f=1:size(cmcs,2)
            r1(f) = cmcs(1,f);
            exps(f) = NM_expectation_from_CMC(cmcs(:,f)');
        end
        weights(i,t,:) = (r1(1:end-1) / 100) / length(r1(1:end-1));
    end
end
    
% Give zero weight to concatnated features
weights = cat(3, weights, nan(size(weights,1), size(weights,2)));

end