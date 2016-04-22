function [issingle] = NM_issingle(data)
% NM_ISSINGLE Check if input data is of SINGLE class type
% 
%   Copyright: Niki Martinel
%   Date: 09/08/2011
%   Return Data: True if input data is of SINGLE class, false otherwise
%   Parameters: 1. any data
%               
%   [ISDOUBLE] = NM_ISSINGLE(DATA) takes input data, DATA, and
%   check if it is of SINGLE class
%
%   DATA should be any kind of data
%   
%   OUT is a TRUE/FALSE if input data is of SINGLE class type or not
%

issingle = strcmpi(class(data), 'single');

end