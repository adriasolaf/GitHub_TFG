function [ value ] = aga_aux_clamp ( value, range )
%aga_aux_clamp Clamp scalar value
%   Limits a scalar value to a closed interval.
%
% Inputs:
%   value: scalar value
%   range: interval [min,max]
%
% Outputs:
%   value: clamped scalar value
%
% Example:
%   [ value ] = aga_aux_clamp ( value, range );
%
% References:
%   [-]
%
%May 2026

    value = max(range(1), min(range(2), value));
end
