function [h] = NM_weightedHistc( inputData, binEdges, weights, excludeRanges, normalize)
%NM_WEIGHTEDHISTC Compute the weighted histogram of input data
%
%   Copyright: Niki Martinel
%   Date: 01/13/2012
%   Return Data: bin elements inside each bin
%   Parameters: 1. input image
%               2. bin edges
%               3. weights
%               4. exclude values in specified ranges
%
%   [H] = NM_WEIGHTEDHISTC(INPUTDATA, BINEDGES, WEIGHTS, EXCLUDERANGE) 
%   takes an input data, INPUT DATA, and the edges of each bin into which put 
%   values, BINEDGES, and produces an ouput value BINELEMENTS 
%   that count the number of elements that are in each bin.
%   Each element of input data is weighted by WEIGHTS matrix before values
%   are put into the histogram
%
%   INPUTDATA, BINEDGES, and EXCLUDE RANGES should be equal to NM_HISTC
%   function
%
%   BINIEDGES should be a vector (which must contain monotonically
%   nondecreasing values)
%
%   EXCLUDERANGES should be a matrix of size Nx2 which rows elements define
%                 ranges to exclude from histogram computation  
%   
%   BINELEMENTS is the number of elements per bin
%

%% Check input parameters
if nargin < 3
    error('nm_weightedhistc:argChk', 'Wrong number of input parameters (<3)');
elseif nargin == 3
    excludeRanges = [];
end

if ~isequal(size(inputData), size(weights))
    error('INPUTDATA and WEIGHTS must be vectors/matrix of the same size');
end

% Exclude ranges from histogram computation
dataForHistogram = inputData(:);
[rows, cols] = size(excludeRanges);
for i=1:rows
    dataForHistogram(excludeRanges(i,1)<=dataForHistogram & dataForHistogram<=excludeRanges(i,2)) = -inf;
end
inputData = dataForHistogram;

% Mex version of weighted histogram
h = single( whistcY(inputData(:), weights(:), binEdges) );

% Do not leave the last bin alone
h(end-1) = h(end-1) + h(end);
h(end) = [];

% Normalize histogram to sum up to 1
if nargin == 5 && normalize
    h = h/norm(h,1);
end

% Due to normalization (if the weights are all 0s) we may end up with NaN
% entries
h(isnan(h)) = 0;

end

