function [tests, finalDists, cmc, featDists, cmcFeat] = test_dist(model, X, testSet, dataset, pars, useLearnedWeights, loadIfExisting)
% Author:    Niki Martinel
% Date:      2014/04/28 09:00:27
% Revision:  0.1
% Copyright: Niki Martinel, 2014

if loadIfExisting
    fprintf('Evaluating models...');
    tTestModel = tic;
end

% Try to load data
testsFile = fullfile(pars.settings.dataFolder, [pars.settings.outputPrefix '_tests.mat']);
if exist(testsFile, 'file') && loadIfExisting
    load(testsFile);
else

    weights = pars.classifier.weights;
    coeffs = pars.classifier.coefficients;
    tests = testSet;
    featDists = cell(size(pars.settings.testCams,1), pars.settings.numTests);
    finalDists = cell(size(pars.settings.testCams,1), pars.settings.numTests);
    cmc = cell(size(pars.settings.testCams,1), pars.settings.numTests);
    cmcFeat = cell(size(cmc));

    % Loop through all camera pairs
    for i=1:size(pars.settings.testCams,1)

        %idxa = find(dataset.cam==pars.settings.testCams(i,1));
        %idxb = find(dataset.cam==pars.settings.testCams(i,2));

        for t=1:pars.settings.numTests

            %s = model(i,t).kissme;
            %probe = unique(testSet(i,t).ID(:,1), 'stable');
            %gallery = unique(testSet(i,t).ID(:,2), 'stable');

            % -----------------------------------------------------------------
            % test on second half
            XX = X;
            if pars.classifier.learnSeparate
                for f=1:length(coeffs)
                    if coeffs(f) > 0 && ~isempty(XX{f})
                        XX{f} = XX{f}(1:coeffs(f),:);
                    else
                        XX{f} = [];
                    end
                end
            end

            names = fieldnames(model(i,t));
            idx = strcmpi(names, 'weights');
            names(idx) = [];
            for nameCounter=1:length(names)       

                % Get weights
                if useLearnedWeights
                    w = model(i,t).weights;
                    % Compute scores only
                    [cmc{i,t}, finalDists{i,t}, cmcFeat{i,t}, featDists{i,t}, ...
                        score, scoreFeat] = my_calcMCMC(model(i,t).(names{nameCounter}), XX, testSet(i,t), w, ...
                        pars.classifier.fusionFunction, true, pars.classifier.multipleShotFusionFunc);

                else
                    w = weights(i,:);

                    % Compute scores and CMCs
                    [cmc{i,t}, finalDists{i,t}, cmcFeat{i,t}, featDists{i,t}, ...
                        score, scoreFeat] = my_calcMCMC(model(i,t).(names{nameCounter}), XX, testSet(i,t), w, ...
                        pars.classifier.fusionFunction, false, pars.classifier.multipleShotFusionFunc);

                end
                
                % Get score
                reshapedScores =  reshape(finalDists{i,t}', length(testSet(i,t).label), []);
                if numel(reshapedScores) == size(tests(i,t).ID,1)
                    tests(i,t).score = score;

                    % Individual features score
                    tests(i,t).scoreFeat = zeros(size(testSet(i,t).ID,1),size(featDists{i,t},3));
                    for f=1:size(featDists{i,t},3)
                        tests(i,t).scoreFeat(:,f) = scoreFeat(:,f);
                    end
                else
                    tests(i,t).score = nan(size(tests(i,t).ID,1),1);

                    for f=1:size(featDists{i,t},3)
                        tests(i,t).scoreFeat(:,f) = nan(size(tests(i,t).ID,1),1);
                    end
                end

            end
        end
    end
    
    % Save data
    try
        if loadIfExisting
            save(testsFile, 'tests', 'finalDists', 'cmc', 'featDists', 'cmcFeat');
        end
    catch ME
        warning('test_dist:saveTests', 'Unable to save tests on file %s.', testsFile)
    end
end

if loadIfExisting
    fprintf('done in %.2f(s)\n', toc(tTestModel));
end

end
