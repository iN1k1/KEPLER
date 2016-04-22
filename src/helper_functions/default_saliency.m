function [saliency] = default_saliency(imageSize, muXY, sigmaXY)
%DEFAULT_SALIENCY Summary of this function goes here
% 
% [OUTPUTARGS] = DEFAULT_SALIENCY(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% Author:    Niki Martinel
% Date:      2014/05/19 12:45:12
% Revision:  0.1
% Copyright: Niki Martinel, 2014

saliency = zeros(imageSize);
[X,Y] = meshgrid(1:imageSize(2), 1:imageSize(1));
XX = ((X-muXY(1)).^2/(sigmaXY(1)^2));
YY = ((Y-muXY(2)).^2/(sigmaXY(2)^2));
saliency = exp( -( XX + YY) );
saliency = repmat(saliency/max(saliency(:)), [1, 1, 3]);

end
