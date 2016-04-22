function [ dataset ] = load_dataset( pars )
%LOAD_DATASET Summary of this function goes here
%   Detailed explanation goes here

t=tic;
fprintf('Load %s dataset...', pars.dataset.name);
namefile = [pars.dataset.name, '.mat'];
load(namefile);
fprintf('done in %.2f(s)\n', toc(t));

end

