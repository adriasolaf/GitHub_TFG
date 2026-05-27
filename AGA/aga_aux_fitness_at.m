function [ fit ] = aga_aux_fitness_at ( fitIn, idx )
%aga_aux_fitness_at Fitness vector safe read
%   Reads a fitness value and returns zero when the requested index is not
%   available.
%
% Inputs:
%   fitIn: input fitness vector
%   idx: requested index
%
% Outputs:
%   fit: selected fitness value
%
% Example:
%   [ fit ] = aga_aux_fitness_at ( fitIn, idx );
%
% References:
%   [-]
%
%May 2026

    if numel(fitIn) >= idx
        fit = fitIn(idx);
    else
        fit = 0;
    end
end
