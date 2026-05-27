function [ decoded, data ] = aga_fun_decode ( ind, ranges )
%aga_fun_decode Decode resonant AGA chromosome
%   Clamps a chromosome to the configured ranges and converts angular genes
%   from degrees to radians for solver calls.
%
% Inputs:
%   ind: AGA individual [vinfo_km_s,theta_deg,phi_deg,nu_deg,lp]
%   ranges: chromosome bounds, one [min,max] row per gene
%
% Outputs:
%   decoded: clamped chromosome [vinfo_km_s,theta_deg,phi_deg,nu_deg,lp]
%   data: decoded chromosome data structure
%
% Example:
%   [ decoded, data ] = aga_fun_decode ( ind, ranges );
%
% References:
%   [-]
%
%May 2026

    if nargin < 2 || isempty(ranges)
        ranges = aga_aux_default_ranges();
    end

    % Force a row vector so every operator sees the same shape.
    decoded = ind(:)';

    % Reject malformed chromosomes early and return a valid-shaped payload.
    if numel(decoded) ~= 5 || any(~isfinite(decoded))
        decoded = [NaN NaN NaN NaN NaN];
        data = aga_aux_empty_data(false);
        return;
    end

    % Clamp continuous genes and force lp to a binary flag.
    decoded(1) = aga_aux_clamp(decoded(1), ranges(1, :));
    decoded(2) = aga_aux_clamp(decoded(2), ranges(2, :));
    decoded(3) = aga_aux_clamp(decoded(3), ranges(3, :));
    decoded(4) = aga_aux_clamp(decoded(4), ranges(4, :));
    decoded(5) = double(round(decoded(5)) ~= 0);

    % Convert angular genes to radians for the trajectory solver.
    data = aga_aux_empty_data(true);
    data.vinfo = decoded(1);
    data.theta = deg2rad(decoded(2));
    data.phi = deg2rad(decoded(3));
    data.nu = deg2rad(decoded(4));
    data.lp = decoded(5);
    data.vinf_out = data.vinfo * [cos(data.phi)*cos(data.theta), cos(data.phi)*sin(data.theta), sin(data.phi)];
end
