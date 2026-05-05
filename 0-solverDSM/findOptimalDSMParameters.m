function [dV_best, nu_best, revs_best, lw_best, lp_best, r_m, v_m_minus, v_m_plus, vinf_in, va] = findOptimalDSMParameters(vinf_out, body, jd0, T_p, N, M, apsis_flag, mu_sun, search_nu_handle, full_output_flag)
%   Finds the optimal DSM for a given outgoing v-infinity vector by
%   searching over all revolution splits (revs_before) and Lambert branches
%   (lw, lp). For each combination, the 1D search over nu_DSM is performed 
%   by an external function handle.
%
% Inputs:
%   vinf_out: outgoing v-infinity vector at departure [km/s]
%   body: identifier for the planetary body
%   jd0: departure epoch
%   T_p: orbital period of the resonant body [s]
%   N: number of spacecraft revolutions
%   M: number of body revolutions
%   apsis_flag: DSM reference apsis (1 = apoapsis, 0 = periapsis)
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   search_nu_handle: handle to a 1D search function with:
%                       x_opt = search_nu_handle(@(x) f(x), a, b)
%   full_output_flag: true = recompute full state vectors for the best
%                      solution, false = return only dV_best
%
% Outputs:
%   dV_best: minimum DSM Delta-V found across all splits/branches [km/s]
%   nu_best: optimal DSM true anomaly change [rad]
%   revs_best: optimal number of revolutions before DSM
%   lw_best: optimal Lambert long-way flag
%   lp_best: optimal Lambert long-period flag
%   r_m: DSM position vector [km]
%   v_m_minus: velocity before DSM [km/s]
%   v_m_plus: velocity after DSM [km/s]
%   vinf_in: arrival v-infinity vector [km/s]
%   va: arrival heliocentric velocity [km/s]
%
% References:
%   [1] PhD thesis, David de la Torre Sangra
%
% See also:
%   computeSingleDSMTransferCost, scan_and_refine_1d, evaluateTotalResonantCost
%
% Adria Sola Foixench
% April 2026

    % Initialize outputs
    dV_best = Inf;
    nu_best = NaN;
    revs_best = NaN;
    lw_best = NaN;
    lp_best = NaN;
    r_m = [NaN NaN NaN];
    v_m_minus = [NaN NaN NaN];
    v_m_plus = [NaN NaN NaN];
    vinf_in = [NaN NaN NaN];
    va = [NaN NaN NaN];

    % Search bounds for nu_DSM
    offset = 1e-3;
    nu_min = 0;
    nu_max = pi - offset;

    sec2days = 1 / 86400;
    tof_total_s = M * T_p;
    jd_f = jd0 + (tof_total_s * sec2days);

    % Position of the planets at departure and arrival
    [r_p0, v_p0] = GetBodyICF(body, jd0, mu_sun, 0);
    [r_pf, v_pf] = GetBodyICF(body, jd_f, mu_sun, 0);
    
    planets_state.r_pf = r_pf;
    planets_state.v_pf = v_pf;
    planets_state.jd_f = jd_f;
    planets_state.tof_total_s = tof_total_s;

    % Initial spacecraft orbit
    v_sc0 = v_p0 + vinf_out;
    [sma, ecc, inc, nu0, argp, raan] = ICF2KEP_O(r_p0, v_sc0, mu_sun);
    
    % If the initial orbit is hyperbolic, skip the search
    if ecc >= 1.0
        return; 
    end
    
    % Mean motion properties for computeDSMStateAndTiming
    n_motion = sqrt(mu_sun / (sma^3));
    E0 = 2.0 * atan(sqrt((1.0 - ecc) / (1.0 + ecc)) * tan(nu0 / 2.0));
    M0 = E0 - ecc * sin(E0);

    % Store everything in the orb_init struct
    orb_init.sma = sma; 
    orb_init.ecc = ecc; 
    orb_init.inc = inc;
    orb_init.argp = argp; 
    orb_init.raan = raan;
    orb_init.n_motion = n_motion; 
    orb_init.M0 = M0;
    orb_init.v_sc0 = v_sc0;

    % Loop over all revolution splits
    % Start from the highest revs_before so the Lambert arc has the fewest
    % multi-revolutions
    for revs_try = (N-1):-1:0

        % Number of Lambert multi-revolutions for this split
        mr_lambert_try = N - revs_try - 1;

        % Long-period flag: only relevant for multi-rev Lambert
        if mr_lambert_try > 0
            lp_max = 1;
        else
            lp_max = 0;
        end

        % Loop over all Lambert branches
        for lw_try = 0:1
            for lp_try = 0:lp_max

                % Build function for this (revs, lw, lp)
                f_nu = @(nu) computeSingleDSMTransferCost(nu, revs_try, lw_try, lp_try, N, apsis_flag, mu_sun, false, orb_init, planets_state);

                % Delegate the 1D search over nu to the external method
                nu_opt = search_nu_handle(f_nu, nu_min, nu_max);

                % Skip if search returned no valid solution
                if isnan(nu_opt)
                    continue;
                end

                % Rapid evaluation
                dV_cand = computeSingleDSMTransferCost(nu_opt, revs_try, lw_try, lp_try,N, apsis_flag, mu_sun, false, orb_init, planets_state);

                % Update global best across all splits and branches
                if dV_cand < dV_best
                    dV_best = dV_cand;
                    nu_best = nu_opt;
                    revs_best = revs_try;
                    lw_best = lw_try;
                    lp_best = lp_try;
                end

            end
        end
    end

    % Recompute best solution with all state vectors
    if full_output_flag && ~isinf(dV_best)
        [dV_best, r_m, v_m_minus, v_m_plus, vinf_in, va] = computeSingleDSMTransferCost(nu_best, revs_best, lw_best, lp_best, N, apsis_flag, mu_sun, true, orb_init, planets_state);
    end

end
