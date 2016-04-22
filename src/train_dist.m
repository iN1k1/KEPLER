function [model] = train_dist(X, dataset, trainSet, pars, optimizeWeights, loadIfExisting)

if loadIfExisting
    fprintf('Training models...');
    tTrainModel = tic;
end

%pars.classifier.coefficients = pca_coeffs;
params.N = dataset.peopleCount;
params.pmetric = 0;
params.lambda = 1;
learn_algs = {...
    LearnAlgoKISSME(params), ...
    };

% Try to load data
modelFile = fullfile(pars.settings.dataFolder, [pars.settings.outputPrefix '_model.mat']);
if exist(modelFile, 'file') && loadIfExisting
    load(modelFile);
else

    coeffs = pars.classifier.coefficients;
    inputX = X;

    % Loop through all camera pairs
    trainTime = tic;
    for i=1:size(pars.settings.testCams,1)

        for t=1:pars.settings.numTests

            trainA = trainSet(i,t).index(:,1)';
            trainB = trainSet(i,t).index(:,2)';
            labels = trainSet(i,t).label';

            % -----------------------------------------------------------------
            % Train Metrics
            for aC=1:length(learn_algs)
                cHandle = learn_algs{aC};
                baseStruct = learnPairwise(cHandle,X{1}(1:2,:), trainA, trainB, labels);
                s = repmat(baseStruct, 1, size(coeffs,2));

                if pars.classifier.learnSeparate
                    %par
                    for f=1:size(coeffs,2)
                        if ~isnan(coeffs(i,f)) && coeffs(i,f) > 0 && ~isempty(X{f})
                            s(f) = learnPairwise(cHandle,X{f}(1:coeffs(i,f),:), trainA, trainB, labels);
                        else
                            if coeffs(i,f) == 0
                                s(f) = learnPairwise(cHandle,X{f}, trainA, trainB, labels);
                            end
                        end
                    end
                else
                    X = inputX;
                    Xtmp = [];
                    for f=1:size(coeffs,2)
                        if ~isnan(coeffs(i,f)) && coeffs(i,f) > 0
                           Xtmp = [Xtmp; X{f}(1:coeffs(i,f),:)];
                        end
                    end
                    X = Xtmp;
                    clear Xtmp;

                    s = learnPairwise(cHandle,X,trainSet(i,t).index(:,1)',  trainSet(i,t).index(:,2)', trainSet(i,t).label');            
                end
            end
            if ~isempty(fieldnames(s))
                model(i,t).(cHandle.type) = s;
            end
        end
    end
    trainTime = toc(trainTime);

    % -----------------------------------------------------------------
    % Learn optimal fusion weights
    if optimizeWeights
        optimWeightsTime = tic;
        weights = optimize_weights(inputX, dataset, pars, model, trainSet); 
        optimWeightsTime = toc(optimWeightsTime);
        for i=1:size(pars.settings.testCams,1)
            for t=1:pars.settings.numTests
                model(i,t).weights = squeeze(weights(i,t,:));
            end
        end
    end
    
    % Save data
    try
        if loadIfExisting
            save(modelFile, 'model');
        end
    catch ME
        warning('train_dist:saveModel', 'Unable to save models on file %s.', modelFile)
    end
end

if loadIfExisting
    fprintf('done in %.2f(s)\n', toc(tTrainModel));
end

end