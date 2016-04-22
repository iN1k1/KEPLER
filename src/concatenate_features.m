function concatenatedFeatures = concatenate_features(features, saliency, pars)
% Author:    Niki Martinel
% Date:      2014/02/21 11:59:14
% Revision:  0.1
% Copyright: Niki Martinel, 2014

fprintf('Concatenating feature vectors...');
t = tic;

% Try to load data
featuresFile = fullfile(pars.settings.dataFolder, [pars.settings.outputPrefix '_concatenated_features.mat']);
if exist(featuresFile, 'file')
    load(featuresFile);
else
    
    fnames = fieldnames(features);
    concatenatedFeatures = cell(length(fnames),1);
    for k=1:length(fnames)
        feat = vertcat(features.(fnames{k}));
        if iscell(feat{1})
            hh1 = cellfun(@(x)(vertcat(x{:})), feat, 'UniformOutput', false);
        else
            hh1 = cellfun(@(x)(vertcat(x)), feat, 'UniformOutput', false);
        end
        if pars.features.(fnames{k}).applySqrt
            concatenatedFeatures{k} = sqrt( cell2mat(hh1')' );
        else
            concatenatedFeatures{k} = cell2mat(hh1')';
        end
    end
    
    % All concat
    concatenatedFeatures(end+1,1) = {[]};
    
    % Save data
    try
        save(featuresFile, 'concatenatedFeatures');
    catch ME
        warning('nm_reid_main:saveConcatenatedFeatures', 'Unable to save concatenated features data on file %s.', featuresFile)
    end
    
end

fprintf('done in %.2f(s)\n', toc(t));
end
