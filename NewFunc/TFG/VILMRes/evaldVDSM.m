function dV_out = evaldVDSM(change_angle, orb, aux, revs_try, mr_lambert_try, lp_max)
%   Evaluate the DSM Delta-V for a specific anomaly.
%   Calculates the Delta-V magnitude for a DSM at a specified true 
%   anomaly by evaluating the associated Lambert problem.
%
% Inputs:
%   change_angle: scanned angle used to determine the DSM true anomaly [rad]
%   orb: structure containing the initial heliocentric Keplerian elements 
%   aux: structure containing auxiliary simulation parameters 
%   revs_try: number of full spacecraft revolutions before the DSM point
%   mr_lambert_try: number of full revolutions for the Lambert transfer arc
%   lp_max: flag to evaluate long-period Lambert branches (0 = short only, 1 = long included)
%
% Outputs:
%   dV_out: computed Delta-V magnitude for the DSM [km/s].
%
% References:
%   [-] n/a
%
% See also:
%   checkValidNu, find_best_lambert, scan_nu_pos, ResonantVILM2
%
% Adrià Solà Foixench
% April 2026

    dV_out = inf;

    % Check geometric validity and get DSM state
    [current_valid, r_m_e, v_m_minus_e, tof_e] = checkValidNu(change_angle, orb, aux, revs_try, mr_lambert_try);

    if ~current_valid
        return;
    end

    % Try all Lambert branches (long-period & way)
    [dV_out, ~, ~, ~] = find_best_lambert(r_m_e, aux.r_pf, tof_e, aux.mu_sun, v_m_minus_e, zeros(1,3), mr_lambert_try, lp_max);

end