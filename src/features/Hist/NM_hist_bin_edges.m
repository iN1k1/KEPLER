function [binEdges] = NM_hist_bin_edges(bins, imageColorSpace)
%NM_HIST_BIN_EDGES Summary of this function goes here
% 
% [OUTPUTARGS] = NM_HIST_BIN_EDGES(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% Author:    Niki Martinel
% Date:      2013/08/19 14:08:48
% Revision:  0.1
% Copyright: Niki Martinel, 2013

% If Lab color space is selected color histogram bins 
% should be set to a range that include negative values. 
% According to the Lab color space definition image can have positive and
% negative values for a* and b*.
if strcmpi(imageColorSpace, 'Lab') == 1
    binEdges{1} = 0:(100/bins(1)):100;
    binEdges{2} = -128:(256/bins(2)):128;
    binEdges{3} = -128:(256/bins(3)):128;
%elseif strcmpi(imageColorSpace, 'LCH') == 1
%    binEdges{1} = 0:(100/bins(1)):100;
%    binEdges{2} = 0:(100/bins(2)):100;
%    binEdges{3} = 0:(360/bins(3)):360;
elseif any(strcmpi(imageColorSpace, {'HSV', 'HSL', 'HSI'})) == 1
    binEdges{1} = 0:(360/bins(1)):360;
    binEdges{2} = 0:(1/bins(2)):1;
    binEdges{3} = 0:(1/bins(3)):1;
elseif strcmpi(imageColorSpace, 'YCbCr') == 1
    binEdges{1} = 16:(219/bins(1)):235;
    binEdges{2} = 16:(224/bins(2)):240;
    binEdges{3} = 16:(224/bins(3)):240;
elseif strcmpi(imageColorSpace, 'Luv') == 1
    binEdges{1} = 0:(100/bins(1)):100;
    binEdges{2} = -134:(256/bins(2)):220;
    binEdges{3} = -140:(256/bins(3)):122;
elseif strcmpi(imageColorSpace, 'YPbPr') == 1
    binEdges{1} = 0:(1/bins(1)):1;
    binEdges{2} = -0.5:(1/bins(2)):0.5;
    binEdges{3} = -0.5:(1/bins(3)):0.5;
elseif strcmpi(imageColorSpace, 'YUV') == 1
    binEdges{1} = 0:(1/bins(1)):1;
    binEdges{2} = -0.5:(1/bins(2)):0.5;
    binEdges{3} = -0.5:(1/bins(3)):0.5;  
%elseif strcmpi(imageColorSpace, 'XYZ') == 1
%    binEdges{1} = 0:(95/bins(1)):95;
%    binEdges{2} = 0:(100/bins(2)):100;
%    binEdges{3} = 0:(110/bins(3)):110;
elseif strcmpi(imageColorSpace, 'opponentRGB') == 1
    binEdges{1} = -1/sqrt(2):(sqrt(2)/bins(1)):1/sqrt(2);
    binEdges{2} = -2/sqrt(6):(4/sqrt(6))/bins(2):2/sqrt(6);
    binEdges{3} = 0:(3/sqrt(3))/bins(3):3/sqrt(3);  
else
    binEdges{1} = 0:(1/bins(1)):1;
    binEdges{2} = 0:(1/bins(2)):1;
    binEdges{3} = 0:(1/bins(3)):1;
end

end
