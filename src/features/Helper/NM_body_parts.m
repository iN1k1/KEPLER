function [imageUpperTorso, imageLowerTorso, imageUpperLegs, imageLowerLegs, ...
    maskUpperTorso, maskLowerTorso, maskUpperLegs, maskLowerLegs] = NM_body_parts(imageRGB, mask, pars, convertToColorSpace)
%IMAGE_PARTS Summary of this function goes here
% 
% [OUTPUTARGS] = IMAGE_PARTS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% Author:    Niki Martinel
% Date:      2013/09/27 11:40:18
% Revision:  0.1
% Copyright: Niki Martinel, 2013

if ~isfield(pars, 'detectTorsoAndLegs')
    pars.detectTorsoAndLegs = false;
end
if ~isfield(pars, 'splitEachBodyPart')
    pars.splitEachBodyPart = false;
end

%------------------------------------------------------------------
% Divide upper and lower body part => then to lower/upper torso...
%------------------------------------------------------------------
if pars.detectTorsoAndLegs
    [torso, legs, ~] = NM_div3parts( imageRGB, mask );  
else
    torso = [1 1 size(imageRGB, 2) size(imageRGB,1)*0.5];
    legs = [1 (size(imageRGB,1)*0.5)+1 size(imageRGB, 2) size(imageRGB,1)];
end

if pars.splitEachBodyPart
    splitPoint = torso(2) + ((torso(4)-torso(2))*0.5);
    upperTorso = [torso(1), torso(2), torso(3), splitPoint];
    lowerTorso = [torso(1), splitPoint+1, torso(3), torso(4)];

    splitPoint = legs(2) + ((legs(4)-legs(2))*0.5);
    upperLegs = [legs(1), legs(2), legs(3), splitPoint];
    lowerLegs = [legs(1), splitPoint+1, legs(3), legs(4)];
end

%------------------------------------------------------------------
% Body parts
%------------------------------------------------------------------

% Body and body parts
imageBodyRGB = imageRGB;
imageUpperTorso = [];
imageUpperTorsoRGB = [];
imageLowerTorso = [];
imageLowerTorsoRGB = [];

imageUpperLegs = [];
imageUpperLegsRGB = [];
imageLowerLegs = [];
imageLowerLegsRGB = [];

maskUpperTorso = [];
maskUpperLegs = [];
maskLowerTorso = [];
maskLowerLegs = [];

if pars.detectTorsoAndLegs
    
    if pars.splitEachBodyPart
        imageUpperTorsoRGB = imageRGB(upperTorso(2):upperTorso(4),upperTorso(1):upperTorso(3),:);
        imageLowerTorsoRGB = imageRGB(lowerTorso(2):lowerTorso(4),lowerTorso(1):lowerTorso(3),:);
        imageUpperLegsRGB  = imageRGB(upperLegs(2):upperLegs(4),upperLegs(1):upperLegs(3),:);
        imageLowerLegsRGB = imageRGB(lowerLegs(2):lowerLegs(4),lowerLegs(1):lowerLegs(3),:);
    else
        imageUpperTorsoRGB = imageRGB(torso(2):torso(4),torso(1):torso(3),:);
        imageUpperLegsRGB  = imageRGB(legs(2):legs(4),legs(1):legs(3),:);
        
        maskUpperTorso = mask(torso(2):torso(4),torso(1):torso(3),:);
        maskUpperLegs = mask(legs(2):legs(4),legs(1):legs(3),:);
    end
    
    if isfield(pars, 'torso') && ~isempty(pars.torso) && any(pars.torso)
        imageUpperTorsoRGB = imresize( imageUpperTorsoRGB, pars.torso(1:2));
        if ~isempty(imageLowerTorsoRGB)
            imageLowerTorsoRGB = imresize( imageLowerTorsoRGB, pars.torso(3:4));
        end
    end
    if isfield(pars, 'legs') && ~isempty(pars.legs) && any(pars.legs)
        imageUpperLegsRGB = imresize( imageUpperLegsRGB, pars.legs(1:2));
        if ~isempty(imageLowerLegsRGB)
            imageLowerLegsRGB = imresize( imageLowerLegsRGB, pars.legs(3:4));
        end
    end
    
    % Dense features
    if ~isempty(pars.blockSize) && ~isempty(pars.blockStep) && all(pars.blockStep~=0)
        cellWidth = pars.blockSize(1);
        cellHeight = pars.blockSize(2);    
        imageUpperTorsoRGB = NM_dense_patches(imageUpperTorsoRGB, cellWidth, cellHeight, pars.blockStep(1), pars.blockStep(2));
        if ~isempty(imageLowerTorsoRGB)
            imageLowerTorsoRGB = NM_dense_patches(imageLowerTorsoRGB, cellWidth, cellHeight, pars.blockStep(1), pars.blockStep(2));
        end
        imageUpperLegsRGB  = NM_dense_patches(imageUpperLegsRGB , cellWidth, cellHeight, pars.blockStep(1), pars.blockStep(2));
        if ~isempty(imageLowerLegsRGB)
            imageLowerLegsRGB = NM_dense_patches(imageLowerLegsRGB , cellWidth, cellHeight, pars.blockStep(1), pars.blockStep(2));
        end
    else
        % Need to convert to cell to keep the same object type as for dense
        % patches
        imageUpperTorsoRGB = {imageUpperTorsoRGB};
        if ~isempty(imageLowerTorsoRGB)
            imageLowerTorsoRGB = {imageLowerTorsoRGB};
        end
        
        imageUpperLegsRGB  = {imageUpperLegsRGB};
        if ~isempty(imageLowerLegsRGB)
            imageLowerLegsRGB = {imageLowerLegsRGB};
        end
    end

else
    % Dense features
    if isfield(pars, 'blockSize') && ~isempty(pars.blockSize) && ~isempty(pars.blockStep) && all(pars.blockSize~=0)
        cellWidth = pars.blockSize(1);
        cellHeight = pars.blockSize(2);
        imageUpperTorsoRGB = NM_dense_patches(imageBodyRGB, cellWidth, cellHeight, pars.blockStep(1), pars.blockStep(2));

    else
        % Get whole image
        imageUpperTorsoRGB = {imageBodyRGB};
    end    
end

% Convert the images or patches into the required color space
if convertToColorSpace
    imageUpperTorso = cellfun(@(image)(NM_colorconverter(image, 'RGB', pars.colorSpace)), imageUpperTorsoRGB, 'UniformOutput', false);
    if ~isempty(imageLowerTorsoRGB)
        imageLowerTorso = cellfun(@(image)(NM_colorconverter(image, 'RGB', pars.colorSpace)), imageLowerTorsoRGB, 'UniformOutput', false);
    end

    if ~isempty(imageUpperLegsRGB)
        imageUpperLegs = cellfun(@(image)(NM_colorconverter(image, 'RGB', pars.colorSpace)), imageUpperLegsRGB, 'UniformOutput', false);
    end
    if ~isempty(imageLowerLegsRGB)
        imageLowerLegs = cellfun(@(image)(NM_colorconverter(image, 'RGB', pars.colorSpace)), imageLowerLegsRGB, 'UniformOutput', false);
    end
else
    imageUpperTorso = imageUpperTorsoRGB;
    imageLowerTorso = imageLowerTorsoRGB;
    imageUpperLegs = imageUpperLegsRGB;
    imageLowerLegs = imageLowerLegsRGB;
end

end
