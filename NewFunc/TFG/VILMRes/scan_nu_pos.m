function nu_opt = scan_nu_pos(N_points, N_refine, change_angle_max, orb, aux, revs_try, mr_lambert_try, lp_max_try)
%   Performs a scan to evaluate the DSM Delta-V across a range of 
%   true anomalies, and then refines the global minimum using a ternary search algorithm.
%
% Inputs:
%   N_points: number of evaluation points for the initial scan
%   N_refine: number of iterations for the ternary search refinement
%   change_angle_max: maximum angle limit for the search [rad]
%   orb: structure containing the initial heliocentric Keplerian elements
%   aux: structure containing auxiliary simulation parameters
%   revs_try: number of full spacecraft revolutions before the DSM point
%   mr_lambert_try: number of full spacecraft revolutions in the Lambert arc
%   lp_max_try: maximum long-period parameter for the Lambert solver (0 or 1)
%
% Outputs:
%   nu_opt: optimal scanned angle with the minimum Delta-V [rad]. 
%
% References:
%   [-] n/a
%
% See also:
%   evaldVDSM, ternary_search_nu, ResonantVILM2
%
% Adrià Solà Foixench
% April 2026


    % Evaluate dV at all N_points
    nu_vec = linspace(0, change_angle_max, N_points); % Nu angle vector
    dv_vec = inf * ones(1, N_points); % dV results initialized to inf

    for k = 1:N_points
        dv_vec(k) = evaldVDSM(nu_vec(k), orb, aux, revs_try, mr_lambert_try, lp_max_try);
    end

    % Skip if no valid solution found in this revs_try
    if all(isinf(dv_vec))
        nu_opt = NaN;
        return;
    end

    % Index of minimum dV
    [~, idx_best] = min(dv_vec);

    nu_opt = ternary_search_nu(idx_best, nu_vec, N_points, N_refine, orb, aux, revs_try, mr_lambert_try, lp_max_try);

end