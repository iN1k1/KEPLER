function [results] = main()
% Copyright: Niki Martinel, 2015

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% INIT PAR POOL
pool = gcp('nocreate');
if isempty(pool)
     parpool(10, 'SpmdEnabled', false);
end

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%   INITALIZE PARAMETERS
pars = init_parameters( '001', 'VIPeR', fileparts(which(mfilename)) );

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%   LOAD DATASET
dataset = load_dataset(pars);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%   FEATURE EXTRACTION
[features, saliency] = extract_features( dataset, pars );

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%   CONCATENATE FEATURE VECTORS
concatenatedFeatures = concatenate_features( features, saliency, pars );
clear features saliency;

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% PCA
pcaFeat = reduce_feature_dimensionality(concatenatedFeatures, pars);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% SPLITS
experiment = [];
splitPars = pars;
splitPars.experiment = experiment;
[trainSet, cvSet, testSet] = split_dataset(dataset, splitPars);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% CROSS VALIDATION
pars.classifier.coefficients = cross_validation(dataset, cvSet, pcaFeat, pars);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% TRAIN
model = train_dist(pcaFeat, dataset, trainSet, pars, pars.classifier.optimalWeights, true);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% TEST
tests = test_dist(model, pcaFeat, testSet, dataset, pars, pars.classifier.optimalWeights, true);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% RESULTS
results = compute_results(tests, pars);

end