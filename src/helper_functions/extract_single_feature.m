function [ features, weightPatches, patchesCenters ] = extract_single_feature( imageRGB, mask, weightImg, featType, featPars, pars )
%EXTRACT_SINGLE_FEATURE Summary of this function goes here
%   Detailed explanation goes here

% Get images
featPars.detectTorsoAndLegs = false;
featPars.splitEachBodyPart = false;
featPars.torso = pars.body.torso;
featPars.legs = pars.body.legs;
if isempty(featPars.blockSize)
    featPars.blockSize = pars.body.blockSize;
    featPars.step = pars.body.blockStep;
end
featuresImage = NM_body_parts(imageRGB, mask, featPars, false);
weightPatches = [];
patchesCenters = [];
if ~isempty(weightImg)
    [weightPatches, patchesCenters] = NM_dense_patches(weightImg , featPars.blockSize(1), featPars.blockSize(2), featPars.blockStep(1), featPars.blockStep(2));
end

% Extract features
features = cell(1, length(featuresImage));
for p=1:length(featuresImage)
    if ~isempty(weightPatches)
        featPars.weight = weightPatches{p};
    end
    features{p} = NM_extractFeatures(featuresImage{p}, featType, featPars);
end

end

