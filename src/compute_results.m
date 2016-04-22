function [results] = compute_results(tests, pars )
% Author:    Niki Martinel
% Date:      2013/02/21 10:38:19
% Revision:  0.1
% Copyright: Niki Martinel, 2013

fprintf('Computing results...');
tComputeRes = tic;
    
% Loop through all tests
for i=1:size(tests,1)
    for t=1:size(tests,2)
        [CMC{i,t}, CMCExpectation(i, :, t), ...
            ROC(i,t), nAUCCMC(i,t), nAUCROC(i,t), ...
            queryPersonsIDs{i,t}, matchingIDs{i,t}] ...
            = NM_CMC_ROC(tests(i,t).ID, tests(i,t).score, ...
            'ROC', false, 'isDist', true, 'statsFunHandle', pars.classifier.multipleShotFusionFunc);
    end

    % Average values
    results(i).CMC = mean(vertcat(CMC{i,:}),1);
    results(i).CMCExpectation = mean(CMCExpectation(i,:,:), 3);
    results(i).nAUCCMC = mean(nAUCCMC(i,:)); 
    
    % Best run
    [~, bestRunIdx] = max(nAUCCMC(i,:));
    results(i).CMCBest = CMC{i,bestRunIdx};
    results(i).nAUCCMCBest = nAUCCMC(i,bestRunIdx);
    results(i).CMCExpectationBest = CMCExpectation(i,:,bestRunIdx);
end

fprintf('done in %.2f(s)\n', toc(tComputeRes));

end
