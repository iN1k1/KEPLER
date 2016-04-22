function [imageDataConverted ] = NM_colorconverter( imageData, inputColorSpace, outputColorSpace )
%
%   NM_COLORCONVERTER Convert RGB image into deisred representation
%
%   Copyright: Niki Martinel
%   Date: 08/11/2011
%   Return Data: return converted image
%   Parameters: 1. Image data to convert
%               2/3. Input color space and output color space
%                   - RGB
%                   - YPbPr
%                   - YCbCr
%                   - JPEG-YCbCr
%                   - YUV
%                   - YIQ
%                   - YDbDr
%                   - HSV or HSB (ranges: 0-360/0-1/0-1)
%                   - HSL or HLS (ranges: 0-360/0-1/0-1)
%                   - HSI (ranges: 0-360/0-1/0-1)
%                   - Lab
%                   - Luv
%                   - LCH
%                   - CAT02 LMS
%

% Check arguments
switch nargin
    case 0 
    case 1
    case 2
        error('nm_colorconverter:argChk', 'Invalid number of input arguments (<3)');
    %case 2
    %    % Use MATALB makecform conversion
    %    cform = makecform(inputColorSpace);
    %    imageDataConverted = applycform(imageData, cform); 
end

if any(strcmpi(inputColorSpace, {'gray','grayscale'}))
    imageData = gray2rgb(imageData);
    inputColorSpace = 'RGB';
end

if strcmpi(outputColorSpace, inputColorSpace) == 1
    imageDataConverted = imageData;
    
elseif any(strcmpi(outputColorSpace, {'gray', 'grayscale'}))
    if strcmpi(inputColorSpace, 'RGB') ~= 1
        imageData = NM_colorconverter(imageData, inputColorSpace, 'RGB');
    end
    imageDataConverted = rgb2gray(imageData);
    
elseif strcmpi(outputColorSpace, 'normRGB') == 1
    
    imageDataConverted(:,:,1) = imageData(:,:,1) ./ (sum(imageData,3) + eps);
    imageDataConverted(:,:,2) = imageData(:,:,2) ./ (sum(imageData,3) + eps);
    imageDataConverted(:,:,3) = sum(imageData,3) ./ (3 + eps);
    
elseif strcmpi(inputColorSpace, 'RGB') && strcmpi(outputColorSpace, 'opponentRGB')
    [R,G,B] = NM_image2channels(imageData);
    imageDataConverted(:,:,1) = (R-G)./sqrt(2);
    imageDataConverted(:,:,2) = (R+G-2*B)./sqrt(6);
    imageDataConverted(:,:,3) = (R+G+B)./sqrt(3);
%elseif strcmpi(inputColorSpace, 'Lab') && strcmpi(outputColorSpace, 'RGB')
%    imageDataConverted = lab2rgb(imageData);
else
    % Use super fast colorspace
    % (http://www.mathworks.com/matlabcentral/fileexchange/28790-colorspace-tra
    % nsformations)
    
    % It requires data in double format or UINT8, so, if the inputa image
    % is single, we just convert the data type do double
    if NM_issingle(imageData)
        imageDataConverted = colorspace([outputColorSpace '<-' inputColorSpace], double(imageData));
        imageDataConverted = single(imageDataConverted);
    else
        imageDataConverted = colorspace([outputColorSpace '<-' inputColorSpace], imageData);
    end
end

end

