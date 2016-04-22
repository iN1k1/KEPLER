function [ features, saliencies ] = extract_features( dataset, pars )

fprintf('Extracting features...');
t = tic;

% Try to load data
featuresFile = fullfile(pars.settings.dataFolder, [pars.settings.outputPrefix '_features.mat']);
if exist(featuresFile, 'file')
    load(featuresFile);
else

    % -------------------------------------------------------------------------
    % Available features
    availableFeatures = fieldnames(pars.features);
    availableFeatures = availableFeatures(structfun(@(x)(x.enabled), pars.features));

    % -------------------------------------------------------------------------
    % SPARSE HAAR FEATURES
    initstate = [1 1 pars.features.sparseHaarHSV.blockSize];
    pars.features.sparseHaarHSV.trparams.initstate = initstate;% object position [x y width height]
    pars.features.sparseHaarHSV.clfparams.width = initstate(3);
    pars.features.sparseHaarHSV.clfparams.height= initstate(4);

    pars.features.sparseHaarLab.trparams.initstate = initstate;% object position [x y width height]
    pars.features.sparseHaarLab.clfparams.width = initstate(3);
    pars.features.sparseHaarLab.clfparams.height= initstate(4);

    %compute feature template
    [ftr.px,ftr.py,ftr.pw,ftr.ph,ftr.pwt] = HaarFtr(pars.features.sparseHaarHSV.clfparams,...
        pars.features.sparseHaarHSV.ftrparams,pars.features.sparseHaarHSV.M);
    
    
    % -------------------------------------------------------------------------
    % init features structure
    features = repmat(cell2struct(repmat({[]},length(availableFeatures),1), availableFeatures), dataset.count, 1);
    saliencies = [];
    images = dataset.images;
    masks = dataset.masks;


    % -------------------------------------------------------------------------
    %    Default Saliency
    % -------------------------------------------------------------------------
    defaultSaliency = ones(size(images(:,:,:,1)));
    if pars.saliency.useDefault
        sz = [size(defaultSaliency,1) size(defaultSaliency,2)];
        muXY = sz([2 1])/2;
        sigmaXY = sz([2 1])/2;
        defaultSaliency = default_saliency(sz,muXY,sigmaXY);
    end

    % Loop through all images
    parfor i=1:size(images,4)
        
        % -----------------------------------------------------------------
        %    Main feature extraction part
        %------------------------------------------------------------------
        image = im2double(images(:,:,:,i));
        origMask = masks(:,:,i);
        mask = logical(ones(size(image,1), size(image,2)));

        % -----------------------------------------------------------------
        %    Default Saliency
        %------------------------------------------------------------------
        saliency = defaultSaliency;

        % -----------------------------------------------------------------
        %    Kernelized Graph-Based Visual Saliency
        %------------------------------------------------------------------
        if isempty(pars.saliency.kernelFunFeats)
            [sal, salimg] = NM_saliency(NM_reid_image_preprocessing(image, mask, 'RGB', 'RGB', 'eq', 0), pars.saliency.method);
        else
            
            % Kernels free parameters
            param = makeGBVSParams;
            param.sigma_frac_act = pars.saliency.kernelFreeParFeats;
            param.sigma_frac_norm = pars.saliency.kernelFreeParMaps;
            
            [sal, salimg] = NM_saliency(NM_reid_image_preprocessing(image, mask, 'RGB', 'RGB', 'eq', 0), pars.saliency.method, ...
                pars.saliency.kernelFunFeats, ...
                pars.saliency.kernelFunMaps, ...
                param);
        end
        saliencies(:,:,:,i) = salimg;
        saliency = pars.saliency.fusionFunction(saliency, salimg);
        saliency = saliency ./ max(saliency(:));

        %------------------------------------------------------------------
        % Extract features
        %------------------------------------------------------------------
        for f=1:length(availableFeatures)

            % Preprocess image
            % The image is projected to the desired color space
            imageToProcess = NM_reid_image_preprocessing( image, mask, pars.dataset.imageColorSpace, ...
                pars.features.(availableFeatures{f}).colorSpace, pars.features.(availableFeatures{f}).histOp, -inf, false);

            % Convert the image color range so that it spans [0,1]
            colorRange = NM_color_range(pars.features.(availableFeatures{f}).colorSpace);
            for ch=1:size(imageToProcess,3)
                imageToProcess(:,:,ch) = (imageToProcess(:,:,ch) - colorRange(ch,1)) / (colorRange(ch,2)-colorRange(ch,1));
            end
            
            % Define feature type
            type = get_feature_type(availableFeatures{f});
            featPars = pars.features.(availableFeatures{f});
            featPars.px = ftr.px;
            featPars.py = ftr.py;
            featPars.pw = ftr.pw;
            featPars.ph = ftr.ph;
            featPars.pwt = ftr.pwt;
            
            % Fix color space to RGB so as the considered histogram range
            % is [0,1]
            featPars.colorSpace = 'RGB';

            % Extract features
            [features(i).(availableFeatures{f})] ...
                = extract_single_feature(imageToProcess, origMask, saliency, type, featPars, pars);
        end   
    end
    
    % Save data
    try
        save(featuresFile, 'features', 'saliencies');
    catch ME
        warning('nm_reid_main:saveFeatures', 'Unable to save features data on file %s.', featuresFile)
    end
end

% Features extraction time
fprintf('done in %.2f(s)\n', toc(t));

end

function type = get_feature_type(feature)
type = feature;
if any(strcmpi(feature, {'HSV', 'HSL', 'HSI', 'RGB', 'normRGB', 'LCH', 'Lab', 'Luv', 'YCbCr', 'YPbPr', 'XYZ', 'YUV'}))
    type = 'whist';
elseif strfind(feature, 'colorMean')
    type = 'colorMean';
elseif strfind(feature, 'sparseHaar')
    type = 'sparseHaar';
elseif strfind(feature, 'sift')
    type = 'sift';
end

end
