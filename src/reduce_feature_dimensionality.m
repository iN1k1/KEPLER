function [pcaFeat] = reduce_feature_dimensionality(features, pars)

fprintf('Reducing features dimensionality (PCA)...');
t = tic;

% Try to load data
dimReducedFeaturesFile = fullfile(pars.settings.dataFolder, [pars.settings.outputPrefix '_dim_reduced_concatenated_features.mat']);
if exist(dimReducedFeaturesFile, 'file')
    load(dimReducedFeaturesFile);
else
    
    pcaFeat = cell(1, length(features));
    if strcmpi(pars.classifier.dimRedMethod, 'pca')
        for i=1:length(features)
            if ~isempty(features{i})
                [pcaFeat{i}, mapping] = pca(features{i}, 200);
                pcaFeat{i} = pcaFeat{i}';
            end
        end
    else
        pcaFeat = cellfun(@(x)(x'), features, 'UniformOutput', false)';
        pars.classifier.coefficients = NaN*ones(size(pars.classifier.weights));
        pars.classifier.coefficients(pars.classifier.weights~=0) = 0;
    end

	% Save data
    try
        save(dimReducedFeaturesFile, 'pcaFeat');
    catch ME
        warning('nm_reid_main:saveDimReducedConcatenatedFeaturesAvg', 'Unable to save  concatenated features with reduced dimensionality on file %s.', dimReducedFeaturesFile)
    end	
end

fprintf('done in %.2f(s)\n', toc(t));

end