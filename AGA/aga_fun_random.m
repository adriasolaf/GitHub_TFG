function [ ind ] = aga_fun_random ( ranges )
%aga_fun_random Random resonant AGA chromosome
%   Generates one bounded individual with format
%   [vinfo_km_s,theta_deg,phi_deg,nu_deg,lp].
%
% Inputs:
%   ranges: chromosome bounds, one [min,max] row per gene
%
% Outputs:
%   ind: random AGA individual
%
% Example:
%   [ ind ] = aga_fun_random ( ranges );
%
% References:
%   [-]
%
%May 2026

    if nargin < 1 || isempty(ranges)
        ranges = aga_aux_default_ranges();
    end

    % Draw continuous genes uniformly inside their configured intervals.
    ind = [
        aga_aux_rand_in_range(ranges(1, :)), ...
        aga_aux_rand_in_range(ranges(2, :)), ...
        aga_aux_rand_in_range(ranges(3, :)), ...
        aga_aux_rand_in_range(ranges(4, :)), ...
        double(rand() >= 0.5)
    ];
end
