function nu_opt = ternary_search_nu(idx_best, nu_vec, N_points, N_refine, orb, aux, revs_try, mr_lambert_try, lp_max_try)
%   Find the exact location of the optimal DSM by performing a 
%   ternary search around the minimum Delta-V found during
%   the initial scan.
%
% Inputs:
%   idx_best: index of the minimum Delta-V found in the initial scan array
%   nu_vec: array of scanned angles from the initial scan [rad]
%   N_points: total number of evaluation points used in the first scan
%   N_refine: number of iterations to perform for the ternary search
%   orb: structure containing the initial heliocentric Keplerian elements
%   aux: structure containing auxiliary simulation parameters
%   revs_try: number of full spacecraft revolutions before the DSM point
%   mr_lambert_try: number of full spacecraft revolutions in the Lambert arc
%   lp_max_try: maximum long-period parameter for the Lambert solver (0 or 1)
%
% Outputs:
%   nu_opt: refined angle with the global minimum Delta-V [rad]
%
% References:
%   [1] Ternary Search (cp-algorithms.com)
%
% See also:
%   scan_nu_pos, evaldVDSM, ResonantVILM2
%
% Adrià Solà Foixench
% April 2026

    % Define refinement bounds
    if idx_best > 1
        nu_lo_r = nu_vec(idx_best - 1);
    else
        nu_lo_r = nu_vec(1);
    end

    if idx_best < N_points
        nu_hi_r = nu_vec(idx_best + 1);
    else
        nu_hi_r = nu_vec(N_points);
    end

    % Refine iterations
    for iter = 1:N_refine
        a = nu_lo_r + (nu_hi_r - nu_lo_r) / 3; % First third point
        b = nu_hi_r - (nu_hi_r - nu_lo_r) / 3; % Second third point

        % Evaluate dV at both third points
        da = evaldVDSM(a, orb, aux, revs_try, mr_lambert_try, lp_max_try);
        db = evaldVDSM(b, orb, aux, revs_try, mr_lambert_try, lp_max_try);

        % Narrow the search interval
        if da <= db
            nu_hi_r = b;
        else
            nu_lo_r = a;
        end
    end

    % Optimal nu angle
    nu_opt = 0.5 * (nu_lo_r + nu_hi_r);

end