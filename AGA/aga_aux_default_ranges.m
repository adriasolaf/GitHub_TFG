function [ ranges ] = aga_aux_default_ranges ()
%aga_aux_default_ranges Default resonant AGA ranges
%   Returns default chromosome bounds. Rows are [vinfo,theta,phi,nu,lp].
%   Angle bounds are in degrees.
%
% Inputs:
%   [-]
%
% Outputs:
%   ranges: chromosome bounds, one [min,max] row per gene
%
% Example:
%   [ ranges ] = aga_aux_default_ranges ();
%
% References:
%   [-]
%
%May 2026

    ranges = [
        0, 20
       -180, 180
       -180, 180
       -180, 180
        0, 1
    ];
end
