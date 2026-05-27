function [ data ] = aga_aux_empty_data ( isValid )
%aga_aux_empty_data Empty decoded chromosome data
%   Builds the decoded-chromosome output template.
%
% Inputs:
%   isValid: chromosome validity flag
%
% Outputs:
%   data: decoded chromosome data structure
%
% Example:
%   [ data ] = aga_aux_empty_data ( isValid );
%
% References:
%   [-]
%
%May 2026

    data = struct('isValid', isValid, 'vinfo', NaN, 'theta', NaN, ...
        'phi', NaN, 'nu', NaN, 'lp', NaN, 'vinf_out', [NaN NaN NaN]);
end
