function [pars] = init_parameters(varargin)
% Author:    Niki Martinel
% Date:      2012/10/30 10:38:19
% Revision:  0.1
% Copyright: Niki Martinel, 2012

fprintf('Initialize parameters...');
t = tic;


p = inputParser;
p.addRequired('testID');
p.addRequired('datasetName');
p.addRequired('folder');
p.parse(varargin{:});
opts = p.Results;

% Test ID
pars.settings.testID = opts.testID;
    
% dataset parameters
pars.dataset.name = opts.datasetName;
pars.dataset.imageColorSpace = 'RGB';

% image sizes
pars.dataset.imageWidth         = 64;
pars.dataset.imageHeight        = 128;
pars.dataset.imageMagFactor     = 1;
pars.dataset.useMasks           = false;  
pars.dataset.loadCams           = false;  

%% ========================================================================
%   BODY PARTS
% =========================================================================
pars.body.detectTorsoAndLegs = false;
pars.body.splitEachBodyPart = false;
pars.body.torso = [64 64 0 0];        % rows (h), cols (w)
pars.body.legs  = [64 64 0 0];
pars.body.blockSize = [16 16];
pars.body.blockStep = [8 8];
pars.body.mag   = 1;


%% ========================================================================
%   SALIENCY
% =========================================================================
pars.saliency.method = 'kgbvs';
pars.saliency.enabled = true;
pars.saliency.useDefault = true;
pars.saliency.fusionFunction = @(X,Y)(X+Y);
pars.saliency.kernelFunFeats = [];
pars.saliency.kernelFreeParFeats = [];
pars.saliency.kernelFunMaps = [];
pars.saliency.kernelFreeParMaps = [];
%pars.saliency.usage = 'fw';
% Gaussian
%pars.classifier.kernelFunFeats = @(x,p)(exp( -1 * x / (2 * p^2)));

% Exponential
%pars.classifier.kernelFunFeats = @(x,p)(exp( -1 * sqrt(x) / (2 *  p^2)));

% Laplacian
%pars.classifier.kernelFunFeats = @(x,p)(exp( -1 * sqrt(x) / p));
    
%% ========================================================================
%   CLASSIFIER
% =========================================================================
pars.classifier.method = 'kissme';
pars.classifier.kfold = 5;
pars.classifier.learnSeparate = true;
pars.classifier.weights = [];
pars.classifier.coefficients = [];
pars.classifier.dimRedMethod = 'pca';
pars.classifier.kernelFun = [];
pars.classifier.crossValidateCoefficients = true;
pars.classifier.fusionFunction = @sum;
pars.classifier.optimalWeights = true;
pars.classifier.multipleShotFusionFunc = @mean; %@(d)(prod(repmat(mean(d), n-length(d),1)) * prod(d));

%% ========================================================================
%   FEATURES
% =========================================================================
distance = 'diff';
histOp = 'none';
applySqrt = true;

% NOTE
% BLOCKS ARE [WIDTH HEIGHT]
% STEPS ARE  [RIGHT BOTTOM]

sparseHaarEnabled = true;
sparseHaarBlockSize = [64 16];
sparseHaarStep = [64 8];

colorMeanEnabled = true;
colorMeanBlock = [8 8];
colorMeanStep = [4 4];

histEnabled = true;
histBins = [24 24 24];
histNormalize = true;
histBlock = [64 16];
histStep = [64 8];
histRoundFunction = [];
histNewMaxValue = [];
histDownsamples = 1;

lbpEnabled = true;
lbpBlockSize = [16 16];
lbpStep = [8 8];

siftEnabled = true;
siftBlocSize = [16 16];
siftStep = [8 8];

%----------------------------------------------------------------
% feature parameters
M = 20; % number of all weaker classifiers, i.e,feature pool
minRect = 2;
maxRect = 4;

pars.features.sparseHaarHSV.trparams.init_postrainrad = 1.0;%radical scope of positive samples
pars.features.sparseHaarHSV.ftrparams.minNumRect = minRect;
pars.features.sparseHaarHSV.ftrparams.maxNumRect = maxRect;
pars.features.sparseHaarHSV.M = M;
pars.features.sparseHaarHSV.colorSpace = 'HSV';
pars.features.sparseHaarHSV.histOp = histOp;
pars.features.sparseHaarHSV.distance = distance;
pars.features.sparseHaarHSV.enabled = sparseHaarEnabled;
pars.features.sparseHaarHSV.blockSize = sparseHaarBlockSize;
pars.features.sparseHaarHSV.blockStep = sparseHaarStep;
pars.features.sparseHaarHSV.applySqrt = false;

pars.features.sparseHaarLab.trparams.init_postrainrad = 1.0;%radical scope of positive samples
pars.features.sparseHaarLab.ftrparams.minNumRect = minRect;
pars.features.sparseHaarLab.ftrparams.maxNumRect = maxRect;
pars.features.sparseHaarLab.M = M;
pars.features.sparseHaarLab.colorSpace = 'Lab';
pars.features.sparseHaarLab.histOp = histOp;
pars.features.sparseHaarLab.distance = distance;
pars.features.sparseHaarLab.enabled = sparseHaarEnabled;
pars.features.sparseHaarLab.blockSize = sparseHaarBlockSize;
pars.features.sparseHaarLab.blockStep = sparseHaarStep;
pars.features.sparseHaarLab.applySqrt = false;

pars.features.sparseHaarYUV.trparams.init_postrainrad = 1.0;%radical scope of positive samples
pars.features.sparseHaarYUV.ftrparams.minNumRect = minRect;
pars.features.sparseHaarYUV.ftrparams.maxNumRect = maxRect;
pars.features.sparseHaarYUV.M = M;
pars.features.sparseHaarYUV.colorSpace = 'YUV';
pars.features.sparseHaarYUV.histOp = histOp;
pars.features.sparseHaarYUV.distance = distance;
pars.features.sparseHaarYUV.enabled = sparseHaarEnabled;
pars.features.sparseHaarYUV.blockSize = sparseHaarBlockSize;
pars.features.sparseHaarYUV.blockStep = sparseHaarStep;
pars.features.sparseHaarYUV.applySqrt = false;

pars.features.sparseHaarNormRGB.trparams.init_postrainrad = 1.0;%radical scope of positive samples
pars.features.sparseHaarNormRGB.ftrparams.minNumRect = minRect;
pars.features.sparseHaarNormRGB.ftrparams.maxNumRect = maxRect;
pars.features.sparseHaarNormRGB.M = M;
pars.features.sparseHaarNormRGB.colorSpace = 'normRGB';
pars.features.sparseHaarNormRGB.histOp = histOp;
pars.features.sparseHaarNormRGB.distance = distance;
pars.features.sparseHaarNormRGB.enabled = sparseHaarEnabled;
pars.features.sparseHaarNormRGB.blockSize = sparseHaarBlockSize;
pars.features.sparseHaarNormRGB.blockStep = sparseHaarStep;
pars.features.sparseHaarNormRGB.applySqrt = false;

pars.features.sparseHaarRGB.trparams.init_postrainrad = 1.0;%radical scope of positive samples
pars.features.sparseHaarRGB.ftrparams.minNumRect = minRect;
pars.features.sparseHaarRGB.ftrparams.maxNumRect = maxRect;
pars.features.sparseHaarRGB.M = M;
pars.features.sparseHaarRGB.colorSpace = 'RGB';
pars.features.sparseHaarRGB.histOp = histOp;
pars.features.sparseHaarRGB.distance = distance;
pars.features.sparseHaarRGB.enabled = sparseHaarEnabled;
pars.features.sparseHaarRGB.blockSize = sparseHaarBlockSize;
pars.features.sparseHaarRGB.blockStep = sparseHaarStep;
pars.features.sparseHaarRGB.applySqrt = false;


%----------------------------------------------------------------
% Color mean
pars.features.colorMeanHSV.enabled = colorMeanEnabled;
pars.features.colorMeanHSV.colorSpace = 'HSV';
pars.features.colorMeanHSV.histOp = 'none';
pars.features.colorMeanHSV.newMaxColorValue = 40;
pars.features.colorMeanHSV.roundFunc = @floor;
pars.features.colorMeanHSV.blockSize = colorMeanBlock;
pars.features.colorMeanHSV.blockStep = colorMeanStep;
pars.features.colorMeanHSV.distance = distance;
pars.features.colorMeanHSV.applySqrt = applySqrt;

pars.features.colorMeanLab.enabled = colorMeanEnabled;
pars.features.colorMeanLab.colorSpace = 'Lab';
pars.features.colorMeanLab.histOp = 'none';
pars.features.colorMeanLab.newMaxColorValue = 40;
pars.features.colorMeanLab.roundFunc = @floor;
pars.features.colorMeanLab.blockSize = colorMeanBlock;
pars.features.colorMeanLab.blockStep = colorMeanStep;
pars.features.colorMeanLab.distance = distance;
pars.features.colorMeanLab.applySqrt = applySqrt;

pars.features.colorMeanYUV.enabled = colorMeanEnabled;
pars.features.colorMeanYUV.colorSpace = 'YUV';
pars.features.colorMeanYUV.histOp = 'none';
pars.features.colorMeanYUV.newMaxColorValue = 40;
pars.features.colorMeanYUV.roundFunc = @floor;
pars.features.colorMeanYUV.blockSize = colorMeanBlock;
pars.features.colorMeanYUV.blockStep = colorMeanStep;
pars.features.colorMeanYUV.distance = distance;
pars.features.colorMeanYUV.applySqrt = applySqrt;

pars.features.colorMeanNormRGB.enabled = colorMeanEnabled;
pars.features.colorMeanNormRGB.colorSpace = 'normRGB';
pars.features.colorMeanNormRGB.histOp = 'none';
pars.features.colorMeanNormRGB.newMaxColorValue = 40;
pars.features.colorMeanNormRGB.roundFunc = @floor;
pars.features.colorMeanNormRGB.blockSize = colorMeanBlock;
pars.features.colorMeanNormRGB.blockStep = colorMeanStep;
pars.features.colorMeanNormRGB.distance = distance;
pars.features.colorMeanNormRGB.applySqrt = applySqrt;

pars.features.colorMeanRGB.enabled = colorMeanEnabled;
pars.features.colorMeanRGB.colorSpace = 'RGB';
pars.features.colorMeanRGB.histOp = 'none';
pars.features.colorMeanRGB.newMaxColorValue = 40;
pars.features.colorMeanRGB.roundFunc = @floor;
pars.features.colorMeanRGB.blockSize = colorMeanBlock;
pars.features.colorMeanRGB.blockStep = colorMeanStep;
pars.features.colorMeanRGB.distance = distance;
pars.features.colorMeanRGB.applySqrt = applySqrt;

%----------------------------------------------------------------------
% SIFT paramters
pars.features.siftHSV.enabled = siftEnabled;
pars.features.siftHSV.colorSpace = 'HSV';
pars.features.siftHSV.histOp = histOp;
pars.features.siftHSV.blockSize = siftBlocSize;
pars.features.siftHSV.blockStep = siftStep;
pars.features.siftHSV.distance = distance;
pars.features.siftHSV.frames = [pars.features.siftHSV.blockSize/2 8/3 0]';
pars.features.siftHSV.onlyDescriptors = true;
pars.features.siftHSV.applySqrt = applySqrt;

pars.features.siftLab.enabled = siftEnabled;
pars.features.siftLab.colorSpace = 'Lab';
pars.features.siftLab.histOp = histOp;
pars.features.siftLab.blockSize = siftBlocSize;
pars.features.siftLab.blockStep = siftStep;
pars.features.siftLab.distance = distance;
pars.features.siftLab.frames = [pars.features.siftLab.blockSize/2 8/3 0]';
pars.features.siftLab.onlyDescriptors = true;
pars.features.siftLab.applySqrt = applySqrt;

pars.features.siftYUV.enabled = siftEnabled;
pars.features.siftYUV.colorSpace = 'YUV';
pars.features.siftYUV.histOp = histOp;
pars.features.siftYUV.blockSize = siftBlocSize;
pars.features.siftYUV.blockStep = siftStep;
pars.features.siftYUV.distance = distance;
pars.features.siftYUV.frames = [pars.features.siftHSV.blockSize/2 8/3 0]';
pars.features.siftYUV.onlyDescriptors = true;
pars.features.siftYUV.applySqrt = applySqrt;

pars.features.siftNormRGB.enabled = siftEnabled;
pars.features.siftNormRGB.colorSpace = 'normRGB';
pars.features.siftNormRGB.histOp = histOp;
pars.features.siftNormRGB.blockSize = siftBlocSize;
pars.features.siftNormRGB.blockStep = siftStep;
pars.features.siftNormRGB.distance = distance;
pars.features.siftNormRGB.frames = [pars.features.siftNormRGB.blockSize/2 8/3 0]';
pars.features.siftNormRGB.onlyDescriptors = true;
pars.features.siftNormRGB.applySqrt = applySqrt;

pars.features.siftRGB.enabled = siftEnabled;
pars.features.siftRGB.colorSpace = 'RGB';
pars.features.siftRGB.histOp = histOp;
pars.features.siftRGB.blockSize = siftBlocSize;
pars.features.siftRGB.blockStep = siftStep;
pars.features.siftRGB.distance = distance;
pars.features.siftRGB.frames = [pars.features.siftRGB.blockSize/2 8/3 0]';
pars.features.siftRGB.onlyDescriptors = true;
pars.features.siftRGB.applySqrt = applySqrt;

%----------------------------------------------------------------------
% Color Histogram

pars.features.HSV.enabled = histEnabled;
pars.features.HSV.colorSpace = 'HSV';            % Lab, LCH, HSV, HSI, RGB, XYZ
pars.features.HSV.histOp = 'none';                    % Norm, eq, norm-eq, none 
pars.features.HSV.bins = histBins;
pars.features.HSV.weights = [];
pars.features.HSV.excludeRange = [-inf -1];
pars.features.HSV.normalize = histNormalize;
pars.features.HSV.distance = distance;
pars.features.HSV.blockSize = histBlock;
pars.features.HSV.blockStep = histStep;
pars.features.HSV.newMaxValue = histNewMaxValue;
pars.features.HSV.roundFunc = histRoundFunction;
pars.features.HSV.downsamples = histDownsamples;
pars.features.HSV.useChannels = [1 1 1];
pars.features.HSV.applySqrt = applySqrt;

pars.features.Lab.enabled = histEnabled;
pars.features.Lab.colorSpace = 'Lab';            % Lab, LCH, HSV, HSI, RGB, XYZ
pars.features.Lab.histOp = 'none';                    % Norm, eq, norm-eq, none 
pars.features.Lab.bins = histBins;
pars.features.Lab.weights = [];
pars.features.Lab.excludeRange = [-inf -129];
pars.features.Lab.normalize = histNormalize;
pars.features.Lab.distance = distance;
pars.features.Lab.blockSize = histBlock;
pars.features.Lab.blockStep = histStep;
pars.features.Lab.newMaxValue = histNewMaxValue;
pars.features.Lab.roundFunc = histRoundFunction;
pars.features.Lab.downsamples = histDownsamples;
pars.features.Lab.applySqrt = applySqrt;

pars.features.YUV.enabled = histEnabled;
pars.features.YUV.colorSpace = 'YUV';
pars.features.YUV.histOp = 'none';
pars.features.YUV.bins = histBins;
pars.features.YUV.excludeRange = [-inf -100];
pars.features.YUV.normalize = histNormalize;
pars.features.YUV.distance = distance;
pars.features.YUV.blockSize = histBlock;
pars.features.YUV.blockStep = histStep;
pars.features.YUV.newMaxValue = histNewMaxValue;
pars.features.YUV.roundFunc = histRoundFunction;
pars.features.YUV.downsamples = histDownsamples;
pars.features.YUV.applySqrt = applySqrt;

pars.features.normRGB.enabled = histEnabled;
pars.features.normRGB.colorSpace = 'normRGB';
pars.features.normRGB.histOp = 'none';
pars.features.normRGB.bins = histBins;
pars.features.normRGB.weights = [];
pars.features.normRGB.excludeRange = [-inf -1];
pars.features.normRGB.normalize = histNormalize;
pars.features.normRGB.distance = distance;
pars.features.normRGB.blockSize = histBlock;
pars.features.normRGB.blockStep = histStep;
pars.features.normRGB.newMaxValue = histNewMaxValue;
pars.features.normRGB.roundFunc = histRoundFunction;
pars.features.normRGB.downsamples = histDownsamples;
pars.features.normRGB.applySqrt = applySqrt;

pars.features.RGB.enabled = histEnabled;
pars.features.RGB.colorSpace = 'RGB';
pars.features.RGB.histOp = 'none';
pars.features.RGB.bins = histBins;
pars.features.RGB.weights = [];
pars.features.RGB.excludeRange = [-inf -1];
pars.features.RGB.normalize = histNormalize;
pars.features.RGB.distance = distance;
pars.features.RGB.blockSize = histBlock;
pars.features.RGB.blockStep = histStep;
pars.features.RGB.newMaxValue = histNewMaxValue;
pars.features.RGB.roundFunc = histRoundFunction;
pars.features.RGB.downsamples = histDownsamples;
pars.features.RGB.applySqrt = applySqrt;

%----------------------------------------------------------------------
% Local Binary Patterns
pars.features.lbp.enabled = lbpEnabled;
pars.features.lbp.colorSpace = 'grayscale';
pars.features.lbp.histOp = 'eq';
pars.features.lbp.patchSize = [];
pars.features.lbp.step = [];
pars.features.lbp.normalizedHistogram = false;
pars.features.lbp.points = 8;
pars.features.lbp.radius = 1;
pars.features.lbp.mapping = 'u2';
pars.features.lbp.blockSize = lbpBlockSize;
pars.features.lbp.blockStep = lbpStep;
pars.features.lbp.applySqrt = applySqrt;


%% ========================================================================
%   SETTINGS
% =========================================================================
pars.settings.testCams = [1 2]; %[1 2; 1 3; 2 3];
pars.settings.numTests = 10;
pars.settings.kfold = 5;
pars.settings.numPersons = [];
pars.settings.numberSamplesPerPersonTraining = [-1 -1];
pars.settings.numberSamplesPerPersonTesting = [-1 -1];
pars.settings.trainAndTestWithSamePersons = false;
pars.settings.testPeopleIDs = [];
pars.settings.learningSets = [0.5 0.5];
pars.settings.extendNumberOfImagesPerPersonPerCamera = 0;
pars.settings.availableColorSpaces = {'grayscale', 'HSV', 'HSL', 'HSI', 'RGB', 'LCH', 'Lab', 'Luv', 'YCbCr', 'YPbPr', 'XYZ', 'YUV', 'normRGB'};
  
% Output file on which save test data
pars.settings.outputPrefix = [pars.dataset.name, '_Id', pars.settings.testID];
pars.settings.rootFolder = opts.folder;
pars.settings.dataFolder = fullfile(pars.settings.rootFolder, pars.settings.outputPrefix);
pars.settings.resultsFolder = fullfile(pars.settings.rootFolder, 'results');

if ~exist(fullfile(pars.settings.dataFolder), 'file')
    mkdir(fullfile(pars.settings.dataFolder));
end

fprintf('done in %.2f(s)\n', toc(t));

end