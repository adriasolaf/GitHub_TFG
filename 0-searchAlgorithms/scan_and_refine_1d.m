function [x_opt, x_vec, f_vec] = scan_and_refine_1d(f, a, b, N_points, N_refine)
%   1D minimization: first performs a grid scan to locate
%   approximately the global minimum, then refines around it using ternary
%   search.
%
% Inputs:
%   f: function handle @(x)
%   a: lower bound of the search interval
%   b: upper bound of the search interval
%   N_points: number of grid points for the scan
%   N_refine: number of iterations for refinement
%
% Outputs:
%   x_opt: refined estimate of the minimum location
%   x_vec: vector of coarse scan points 
%   f_vec: vector of coarse scan values
%
% References:
%   [-] n/a
%
% See also:
%   grid_scan_1d
%
% Adria Sola Foixench
% April 2026

    % Coarse grid scan over [a, b]
    [~, ~, x_vec, f_vec] = grid_scan_1d(f, a, b, N_points);

    % Check if any valid solution was found
    if all(isinf(f_vec))
        x_opt = NaN;
        return;
    end

    % Find the grid cell containing the minimum
    [~, idx_best] = min(f_vec);

    % Define refinement bounds as the neighboring grid cell
    if idx_best > 1
        ref_lo = x_vec(idx_best - 1);
    else
        ref_lo = x_vec(1);
    end

    if idx_best < N_points
        ref_hi = x_vec(idx_best + 1);
    else
        ref_hi = x_vec(N_points);
    end

    %Refinement within the cell
    x_fine = linspace(ref_lo, ref_hi, N_refine);
    f_fine = arrayfun(f, x_fine);
    [min_val, idx_fine] = min(f_fine);
    if ~isfinite(min_val)
        x_opt = NaN; 
        return;
    end
    x_opt = x_fine(idx_fine);

end
