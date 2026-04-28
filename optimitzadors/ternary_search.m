function x_opt = ternary_search(f, a, b, N_iter)
%   Finds the minimum of a function within [a, b] using the
%   ternary search algorithm. At each iteration, the interval is reduced
%   by a factor of 2/3 by evaluating f at two interior points.
%
% Inputs:
%   f: function handle @(x)
%   a: lower bound of the search interval
%   b: upper bound of the search interval
%   N_iter: number of refinement iterations
%
% Outputs:
%   x_opt: estimated location of the minimum (midpoint of final interval)
%
% References:
%   [1] Ternary Search, cp-algorithms.com
%
% See also:
%   grid_scan_1d, scan_and_refine_1d
%
% Adria Sola Foixench
% April 2026

    for iter = 1:N_iter
        % Two interior points dividing [a, b] into three equal segments
        m1 = a + (b - a) / 3;
        m2 = b - (b - a) / 3;

        % Evaluate function at both points
        f1 = f(m1);
        f2 = f(m2);

        % Discard the third that cannot contain the minimum
        if f1 <= f2
            b = m2;
        else
            a = m1;
        end
    end

    % Return midpoint of the final interval as the best estimate
    x_opt = 0.5 * (a + b);

end
