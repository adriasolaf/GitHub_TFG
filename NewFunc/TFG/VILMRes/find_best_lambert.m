function [dV_opt, vinf_in_opt, va_opt, v_m_plus_opt] = find_best_lambert(r_m_opt, r_pf, tof_opt, mu_sun, v_m_minus_opt, v_pf, mr_lambert_try, lp_max_try)
%   Solves the Lambert boundary value problem over all valid 
%   long/short-way and long/short-period transfer branches to find the 
%   trajectory that minimizes the required DSM Delta-V.
%
% Inputs:
%   r_m_opt: heliocentric position vector at the optimal DSM point [km]
%   r_pf: heliocentric position vector of the arrival planet [km]
%   tof_opt: optimal time of flight for the Lambert transfer arc [s]
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   v_m_minus_opt: spacecraft heliocentric velocity just before the DSM [km/s]
%   v_pf: heliocentric velocity vector of the arrival planet [km/s]
%   mr_lambert_try: number of full spacecraft revolutions for the Lambert solver
%   lp_max_try: maximum long-period parameter (0 = short period only, 1 = check both)
%
% Outputs:
%   dV_opt: minimum calculated Delta-V magnitude for the DSM [km/s]
%   vinf_in_opt: optimal hyperbolic excess velocity vector at arrival [km/s]
%   va_opt: optimal heliocentric arrival velocity of the spacecraft [km/s]
%   v_m_plus_opt: optimal spacecraft heliocentric velocity just after the DSM [km/s]
%
% References:
%   [-] n/a
%
% See also:
%   Lambert, evaldVDSM, ResonantVILM2
%
% Adrià Solà Foixench
% April 2026

    dV_opt = inf;
    vinf_in_opt = [NaN NaN NaN];
    va_opt = [NaN NaN NaN];
    v_m_plus_opt = [NaN NaN NaN];

    for lp = 0:lp_max_try
        for lw = 0:1
            [v_m_plus_tbm, v_scf_opt, flag_tbm, ~] = Lambert(r_m_opt, r_pf, tof_opt, mu_sun, lw, mr_lambert_try, lp);

            if flag_tbm == 0
                dV_tbm = norm(v_m_plus_tbm - v_m_minus_opt); % Best maneuver try Delta-V

                if dV_tbm < dV_opt
                    dV_opt = dV_tbm;
                    vinf_in_opt = v_scf_opt - v_pf;
                    va_opt = v_scf_opt;
                    v_m_plus_opt = v_m_plus_tbm;
                end
            end
        end
    end

end