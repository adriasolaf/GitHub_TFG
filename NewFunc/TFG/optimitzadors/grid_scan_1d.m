function [x_best, f_best, x_vec, f_vec] = grid_scan_1d(f, a, b, N_points)
%   Evaluates a scalar function over a 1D grid and returns the
%   point with the minimum value. Also returns all evaluated points and
%   values for analysis and plotting.
%
% Inputs:
%   f: function handle @(x)
%   a: lower bound of the search interval
%   b: upper bound of the search interval
%   N_points: number of equally spaced evaluation points
%
% Outputs:
%   x_best: point with the minimum function value
%   f_best: minimum function value found
%   x_vec: vector of all evaluated points 
%   f_vec: vector of all function values
%
% References:
%   [-] n/a
%
% See also:
%   ternary_search, scan_and_refine_1d
%
% Adria Sola Foixench
% April 2026

    % Build uniform grid over [a, b]
    x_vec = linspace(a, b, N_points);
    f_vec = inf(1, N_points);

    % Evaluate function at every grid point
    for k = 1:N_points
        f_vec(k) = f(x_vec(k));
    end

    % Find the grid point with minimum value
    [f_best, idx_best] = min(f_vec);
    x_best = x_vec(idx_best);

end
