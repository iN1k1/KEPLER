function [saliency, salmap, timeToCompute] = NM_saliency(image, method, varargin)
[h,w,~] = size(image);
if h<128
    image = imresize(image, [128 NaN]);
end
if w<128
	image = imresize(image, [NaN 128]);
end

t = tic;
saliency = kgbvs(image, varargin{:});
saliency.final_resized_map = repmat(imresize(saliency.master_map_resized, [h w]), [1 1 3]); 
salmap = saliency.final_resized_map;

% Rescale to [0,1]
salmap = salmap - min(salmap(:));
salmap = salmap / max(salmap(:));

timeToCompute = toc(t);
end

