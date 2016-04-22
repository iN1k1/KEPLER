function [image] = NM_reid_image_preprocessing( image, mask, inputColorSpace, ...
                                                outputColorSpace, histOp, ...
                                                maskOuterValue, convertToColorThenPreProcess)
% Author:    Niki Martinel
% Date:      2012/10/10 08:30:47
% Revision:  0.1
% Copyright: Niki Martinel, 2012

if nargin<7
    convertToColorThenPreProcess = true;
end
  
% Convert image
if convertToColorThenPreProcess
    image = NM_colorconverter(image, inputColorSpace, outputColorSpace);
end

if ~convertToColorThenPreProcess
    image = NM_colorconverter(image, inputColorSpace, outputColorSpace);
end

end
