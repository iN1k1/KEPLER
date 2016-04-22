function [expectation] = NM_expectation_from_CMC( cmc )
%#codegen
coder.inline('never')

% CMC is not in percentages
if any(cmc>100)
    cmc = (cmc ./ sum(cmc, 2)) * 100;
end

% Calculate PDF (prob. density function) by applying a discrete
% deriviative. Divide by 100 because the CMC is in percents.
PDF = ([cmc 0] - [0 cmc]) ./ 100;
PDF(end) = []; % Remove the last element


% Calculate the expectation
expectation = sum( [1:length(cmc)] .* PDF );
end