function [x_best, f_best] = grid_scan_3d(f, grid1, grid2, grid3)
%   Evaluates a function over a 3D cartesian grid and returns the
%   point with the minimum value.
% Inputs:
%   f: function handle @([x1, x2, x3])
%   grid1: vector of values for the first dimension
%   grid2: vector of values for the second dimension
%   grid3: vector of values for the third dimension
%
% Outputs:
%   x_best: 3-element vector [x1, x2, x3] with the minimum function value
%   f_best: minimum function value found
%
% References:
%   [-] n/a
%
% See also:
%   grid_scan_1d, scan_and_refine_1d
%
% Adria Sola Foixench
% April 2026

    f_best = inf;
    x_best = [NaN NaN NaN];

    % Sweep over all combinations of the three grids
    for x1 = grid1
        for x2 = grid2
            for x3 = grid3
                x_cand = [x1, x2, x3];
                f_cand = f(x_cand);

                if f_cand < f_best
                    f_best = f_cand;
                    x_best = x_cand;
                end
            end
        end
    end

end
