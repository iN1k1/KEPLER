function [ extractedFeature, times ] = NM_extractFeatures( image, featureType, algorithmParameters)

%% Check input parameters
if nargin ~= 3
   error('nm_extractfeatures:argCheck', 'Wrong number of input arguments (!=3)');  
end

%% Extract features from image
switch featureType
    
    %   -------------------------------------------------------------------
    %       SIFT
    %   -------------------------------------------------------------------
    case 'sift'
        
        %Sample
        %   points = #of sift points to retrieve = 50
        %   levels =  Set the number of levels per octave of the DoG scale
        %   space = 5
        %   displayImage = true/false = if you want to see or not the
        %   original imag
        %   plotFrame = true/false if you want to see or not the frame kypoinsts = true
        %   plotDescriptors = true/false if you want to see or not the sift
        %                   descriptors = true
        %   colorMeanAndHist = true/false if you want to extract the mean 
        %               and the histogram of some points around the sift frame
        %   colorRadius = radius of the circle on which evaluate color mean
        %               or color histogram (have to be integer, otherwise it is rounded to the nearest one)
        %   colorHistBin = edges of histogram bins
        %   colorHistNormalize = normalize histograms to [0,1]
        %   colorHistImageColorSpace = image color space from which
        %                              histogram are extracted
        %   imageForColorExtraction = image on which perform color
        %                             extraction (if not set use original
        %                             one)

        siftTimer1 = tic;
        
        % Input image must be single class and grayscale
        for i=1:size(image,3)
            imageDataToProcess = single(image(:,:,i));

            % Check if "levels" par is defined
            if ~isfield(algorithmParameters, 'levels')
                algorithmParameters.levels = 3;
            end

            [extractedFeature(i).frameKeypoints, extractedFeature(i).descriptors] = vl_sift(imageDataToProcess, 'frames', algorithmParameters.frames);
            
        end
        times = toc(siftTimer1);

        if algorithmParameters.onlyDescriptors
            tmpFeat = [extractedFeature.descriptors];
            clear extractedFeature;
            for i=1:size(tmpFeat,2)
                extractedFeature{i} = double(tmpFeat(:,i));
            end
            clear tmpFeat;
        end
        
    %   -------------------------------------------------------------------
    %       WEIGHTED HISTOGRAM
    %   -------------------------------------------------------------------
    case 'whist'
        if ~isfield(algorithmParameters, 'excludeRange'),   algorithmParameters.excludeRange = []; end
        if ~isfield(algorithmParameters, 'normalize'),      algorithmParameters.normalize = false; end
        if ~isfield(algorithmParameters, 'weight'),         algorithmParameters.weight = ones(size(image)); end
        if ~isfield(algorithmParameters, 'newMaxValue'),    algorithmParameters.newMaxValue = []; end
        if ~isfield(algorithmParameters, 'roundFunc'),      algorithmParameters.roundFunc = []; end
        if ~isfield(algorithmParameters, 'downsamples'),    algorithmParameters.downsamples = 1; end
        if ~isfield(algorithmParameters, 'useChannels'),    algorithmParameters.useChannels = ones(1, size(image,3)); end
        
        % Histogram bin edges
        colorHistBin = NM_hist_bin_edges(algorithmParameters.bins, algorithmParameters.colorSpace);
        
         % Compute WEIGHTED HISTOGRAM
        histTimer1 = tic;
        fidx = 1;
        for i=1:size(image,3)
            if algorithmParameters.useChannels(i) == 1
                for j=1:algorithmParameters.downsamples
                    weight = algorithmParameters.weight(:,:,i);
                    im = image(:,:,i);
                    bins = colorHistBin{i};
                    if j > 1
                        weight = imresize(weight, 1/(2^(j-1)));
                        im = imresize(im, 1/(2^(j-1)));
                        %bins = bins(1:j:end);
                    end

                    extractedFeature{fidx} = NM_weightedHistc(im, ...
                        bins, ...
                        weight, ...
                        algorithmParameters.excludeRange, ...
                        algorithmParameters.normalize);

                    % Scale to new range
                    if ~isempty(algorithmParameters.newMaxValue)
                        extractedFeature{fidx} = algorithmParameters.newMaxValue * (extractedFeature{fidx} / sum(weight(:)));
                    end

                    % Round extracted features if required
                    if ~isempty(algorithmParameters.roundFunc)
                        extractedFeature{fidx} = algorithmParameters.roundFunc(extractedFeature{fidx});
                    end
                    fidx = fidx + 1;
                end
            end
        end
        times(1) = toc(histTimer1);
        
    %   -------------------------------------------------------------------
    %       LBP
    %   -------------------------------------------------------------------   
    case 'lbp'
        % Compute the Local Binary Pattern (LBP) of the given image
        %   Parameters:
        %       patchSize [] = dimension of the patch (patchSize x patchSize (
        %       if set to [] compute the LBP for the whole image
        %       step [] = step size between one grid sampling and the next one
        %       normalizedHistogram [false] = normalize each histogram in the
        %                                         sampling grid
        %       mapping ['u2'] = type of mapping. See LBP.m 
        %                            'u2' = for uniform LBP
        %                            'ri' = for rotation-invariant LBP
        %                            'riu2' = uniform rotation-invariant LBP.
        %       points [8] = The LBP codes are computed using POINTS sampling 
        %                        points on a circle of radius RADIUS
        %       radius [1] = Sampling points circle radius
        %
        %   Example:
        %   im = imread('cameraman.tif');
        %   lbpPars.patchSize = 8;
        %   lbpPars.step = 4;
        %   lbPars.points = 8;
        %   lbpPars.radius = 1;
        %   lbpPars.mapping = 'riu2';
        %   lbp =  NM_extractFeatures(im, 'lbp', lbpPars)
        %   stem(lbp);
        if ~isfield(algorithmParameters, 'patchSize'),          algorithmParameters.patchSize = []; end
        if ~isfield(algorithmParameters, 'step'),               algorithmParameters.step = []; end
        if ~isfield(algorithmParameters, 'normalizedHistogram'),algorithmParameters.normalizedHistogram = false; end 
        if ~isfield(algorithmParameters, 'mapping'),            algorithmParameters.mapping = 'u2'; end
        if ~isfield(algorithmParameters, 'points'),             algorithmParameters.points = 8; end
        if ~isfield(algorithmParameters, 'radius'),             algorithmParameters.radius = 1; end
        
        if algorithmParameters.normalizedHistogram
            histType = 'nh';
        else
            histType = 'h';
        end
        
        % LBP Mapping
        mapping = getmapping(algorithmParameters.points, algorithmParameters.mapping); 
        extractedFeature = [];
        if isempty(algorithmParameters.patchSize) || isempty(algorithmParameters.step)
            if isempty(algorithmParameters.radius)
                extractedFeature = LBP(image);
            else
                extractedFeature = LBP(image, algorithmParameters.radius, algorithmParameters.points, mapping, histType)';
            end
        else
            % Dense LBP
            tmpAlgorithmParameters = algorithmParameters;
            tmpAlgorithmParameters.feature = 'lbp';
            extractedFeature = NM_extractFeatures(image, 'dense', tmpAlgorithmParameters);
        end
        
    %   -------------------------------------------------------------------
    %       Color Mean
    %   -------------------------------------------------------------------   
    case 'colorMean'
        
        % Color range
        colorRange = NM_color_range(algorithmParameters.colorSpace);
        
        % pars
        if ~isfield(algorithmParameters, 'minMaxColorValue'), algorithmParameters.minMaxColorValue = colorRange; end
        if ~isfield(algorithmParameters, 'newMaxColorValue'), algorithmParameters.newMaxColorValue = 1; end
        if ~isfield(algorithmParameters, 'roundFunc'), algorithmParameters.roundFunc = []; end
        
        % For all colors
        extractedFeature = zeros(size(image,3),1);
        for i=1:size(image,3)
            
            % Exclude inf or nan entries
            data = image(:, :, i);
            data(isinf(data) | isnan(data)) = [];
            imgSize = numel(data);
            
            % Sum of all the pixel values in the given image
            ss = sum(sum(data - algorithmParameters.minMaxColorValue(i,1)));
            
            % Maxium value
            maxValue = imgSize * (algorithmParameters.minMaxColorValue(i,2) - algorithmParameters.minMaxColorValue(i,1));
                
            % Compute the average pixel value
            % to get ss/maxValue in [0,1]
            % then scale to the new max color value
            extractedFeature(i) = algorithmParameters.newMaxColorValue * (ss / maxValue);

        end
        
        % Round extractede features if required
        if ~isempty(algorithmParameters.roundFunc)
            extractedFeature = algorithmParameters.roundFunc(extractedFeature);
        end
            
    %   -------------------------------------------------------------------
    %       Sparse Haar
    %   -------------------------------------------------------------------   
    case 'sparseHaar'
        sampleImage.sx = 1;
        sampleImage.sy = 1;
        sampleImage.sw = size(image,2);
        sampleImage.sh = size(image,1);

        %-----------------------------------
        %--------Feature extraction
        for i=1:size(image,3)
            iH = integral(image(:,:,i)*255);%Compute integral image
            extractedFeature{i} = getFtrVal(iH,sampleImage,algorithmParameters);
        end
end


% Sqrt on features
if isfield(algorithmParameters, 'applySqrt') && algorithmParameters.applySqrt
    if iscell(extractedFeature)
        extractedFeature = cellfun(@(feat)(sqrt(feat)), extractedFeature, 'UniformOutput', false);
    else
        extractedFeature = sqrt(extractedFeature);
    end
end

% L1 / L2 normalize the feature
if (isfield(algorithmParameters, 'L1norm') && algorithmParameters.L1norm) ||  (isfield(algorithmParameters, 'L2norm') && algorithmParameters.L2norm)
    
    norm_type = 1;
    if isfield(algorithmParameters, 'L2norm')
        norm_type = 2;
    end
    if iscell(extractedFeature )
        extractedFeature = cellfun(@(feat)(normc_safe(feat,norm_type)), extractedFeature, 'UniformOutput', false);
    else
        extractedFeature = normc_safe(extractedFeature,norm_type);
    end
end

% Warning if NAN/Inf feature
nan_or_inf = [];
if iscell(extractedFeature )
    if iscell(extractedFeature{1})
        nan_or_inf = cell2mat([extractedFeature{:}]);
        nan_or_inf = isnan(nan_or_inf(:)) | isinf(nan_or_inf(:));
    else
        nan_or_inf = cellfun(@(feat)(any(isnan(feat(:)) | isinf(feat(:)))), extractedFeature, 'UniformOutput', false);
        nan_or_inf = [nan_or_inf{:}];
    end    
else
    nan_or_inf = isnan(extractedFeature(:)) | isinf(extractedFeature(:));
end
if any(nan_or_inf(:))
    warning('NM_extractFeatures:nan_inf', sprintf('NaN/Inf %s Feature', featureType));
end

end
