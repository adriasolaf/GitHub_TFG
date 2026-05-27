function [ value ] = aga_aux_rand_in_range ( range )
%aga_aux_rand_in_range Random scalar in range
%   Draws one uniformly distributed random scalar inside [min,max].
%
% Inputs:
%   range: interval [min,max]
%
% Outputs:
%   value: random scalar value
%
% Example:
%   [ value ] = aga_aux_rand_in_range ( range );
%
% References:
%   [-]
%
%May 2026

    value = range(1) + (range(2) - range(1)) * rand();
end
