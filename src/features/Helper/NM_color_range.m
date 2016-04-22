function [range] = NM_color_range(imageColorSpace)
%NM_COLOR_RANGE Summary of this function goes here
% 
% [OUTPUTARGS] = NM_COLOR_RANGE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% Author:    Niki Martinel
% Date:      2014/02/06 14:55:46
% Revision:  0.1
% Copyright: Niki Martinel, 2014

if strcmpi(imageColorSpace, 'Lab') == 1
    range = [0 100;
            -128 128;
            -128 128];
%elseif strcmpi(imageColorSpace, 'LCH') == 1
%    binEdges{1} = 0:(100/bins(1)):100;
%    binEdges{2} = 0:(100/bins(2)):100;
%    binEdges{3} = 0:(360/bins(3)):360;
elseif any(strcmpi(imageColorSpace, {'HSV', 'HSL', 'HSI'})) == 1
    range = [0 360;
             0 1;
             0 1];
elseif strcmpi(imageColorSpace, 'YCbCr') == 1
    range = [16 235;
             16 240;
             16 240];
elseif strcmpi(imageColorSpace, 'Luv') == 1
    range = [0 100;
             -134 220;
             -14 122];
elseif strcmpi(imageColorSpace, 'YPbPr') == 1
    range = [0 1;
             -0.5 0.5;
    	     -0.5 0.5];
elseif strcmpi(imageColorSpace, 'YUV') == 1
    range = [0 1;
             -0.5 0.5;
    	     -0.65 0.65]; 
%elseif strcmpi(imageColorSpace, 'XYZ') == 1
%    range{1} = [0 (95/bins(1)):95;
%    range{2} = [0 (100/bins(2)):100;
%    range{3} = [0 (110/bins(3)):110;
else
    range = [0 1;
             0 1;
             0 1];
end


end
