function [ coeffs ] = cross_validation(dataset, cvSet, X, pars )
%CROSS_VALIDATION Summary of this function goes here
%   Detailed explanation goes here

fprintf('Performing Cross-Validation...');
t = tic;
 
% Try to load data
cvParamsFile = fullfile(pars.settings.dataFolder, [pars.settings.outputPrefix '_cross_val_coeffs.mat']);
if exist(cvParamsFile, 'file')
    load(cvParamsFile);
else
    % Cross validation
    if ~isempty(cvSet)

        % PCA coefficients list
        cvCoeff = 25:55; % We used 10:70 in the paper..
        coeffs = zeros(size(pars.settings.testCams,1),size(X,2));
        %par
        for jj=1:size(X,2)

            if isempty(X{jj})
                coeffs(:,jj) = 0;
            end
            
            % Expectation and rank1 accumulators
            accExpValues = zeros( size(pars.settings.testCams,1), length(cvCoeff) );
            accRank1Values = zeros( size(pars.settings.testCams,1), length(cvCoeff) );

            % run over the kfolds
            for kf=1:pars.classifier.kfold
                trainSet = cvSet.train(:,:,kf);
                testSet = cvSet.test(:,:,kf);

                % Train models
                %tic
                parfor c=1:length(cvCoeff)

                    try
                        % Number of coeffs is larger than the actual number of
                        % features..
                        if cvCoeff(c) > size(X{jj},1)
                            accExpValues(:,c) = nan(size(pars.settings.testCams,1), 1);
                            accRank1Values(:,c) = nan(size(pars.settings.testCams,1), 1);
                            continue;
                        end

                        % Cross validation parameters
                        cvPars = pars;
                        cvPars.classifier.weights = zeros(1,size(X,2));
                        cvPars.classifier.weights(jj) = 1;
                        cvPars.classifier.coefficients = nan(1,size(X,2));
                        cvPars.classifier.coefficients(jj) = cvCoeff(c);

                        % Train and test..
                        % Get expectation and rank 1 values
                        [rank1, expectation] = train_and_test_dist(X, dataset, trainSet, testSet, cvPars);

                        % Accumulate the results
                        accExpValues(:,c) = accExpValues(:,c) + expectation;
                        accRank1Values(:,c) = accRank1Values(:,c) + rank1;
                    catch er
                        disp(er);
                    end
                        
                end
                %toc
            end

            % Turn the accumulated values to averages
            accExpValues = accExpValues / pars.classifier.kfold;
            accRank1Values = accRank1Values / pars.classifier.kfold;

            % Find the best params for the whole outer-fold
            [bestExpValue, bestExpValueIdx] = min(accExpValues,[],2);
            [bestRank1Value, bestRankValueIdx] = max(accRank1Values, [], 2);

            % Optimal coeff
            coeffs(:,jj) = cvCoeff(bestRankValueIdx);
            coeffs(:,jj) = cvCoeff(bestExpValueIdx);

        end
    end

    % Save data
    try
        save(cvParamsFile, 'coeffs');
    catch ME
        warning('nm_reid_main:saveFeatures', 'Unable to save features data on file %s.', featuresFile)
    end
end

% Features extraction time
fprintf('done in %.2f(s)\n', toc(t));


end


function[rank1, cmcExp] = train_and_test_dist(X, dataset, trainSet, testSet, pars)

[model] = train_dist(X, dataset, trainSet, pars, false, false);
[tests, finalDists, CMC, featDists, featCMC] = test_dist(model, X, testSet, dataset, pars, false, false);
rank1 = zeros(size(pars.settings.testCams,1),1);
cmcExp = zeros(size(pars.settings.testCams,1),1);
for i=1:size(pars.settings.testCams)
    cmc = mean(vertcat(CMC{i,:}));
    rank1(i) = cmc(1);
    cmcExp(i) = NM_expectation_from_CMC(cmc);
end
end