function [current_valid, r_m_out, v_m_minus_out, tof_lambert_out, v_sc0] = checkValidNu(nu_DSM, revs_before, N, apsis_flag, mu_sun, orb_init, planets_state)
%   Computes the spacecraft state at a DSM point and determines
%   if the resulting geometry is valid for a Lambert transfer to the
%   arrival planet. This function performs all the Kepler propagation
%   and validity checks without solving Lambert problem.
%
% Inputs:
%   vinf_out: outgoing v-infinity vector at departure [km/s]
%   nu_DSM: scanned angle for DSM true anomaly placement [rad]
%   revs_before: number of full Kepler revolutions before the DSM
%   body: identifier for the planetary body
%   jd0: departure epoch
%   T_p: orbital period of the resonant body [s]
%   N: number of spacecraft revolutions
%   M: number of body revolutions
%   apsis_flag: DSM reference apsis (1 = apoapsis, 0 = periapsis)
%   mu_sun: central body gravitational parameter [km^3/s^2]
%
% Outputs:
%   current_valid: boolean, true if the DSM geometry is valid
%   r_m_out: heliocentric position vector at the DSM point [km]
%   v_m_minus_out: heliocentric velocity vector just before the DSM [km/s]
%   tof_lambert_out: remaining time of flight for the Lambert arc [s]
%   v_sc0: spacecraft heliocentric velocity at departure [km/s]
%
% References:
%   [-] n/a
%
% See also:
%   build_RessOrbDSM, checkLambertValidity
%
% Adria Sola Foixench
% April 2026

    % Initialize outputs
    current_valid = false;
    r_m_out = zeros(1,3);
    v_m_minus_out = zeros(1,3);
    tof_lambert_out = 0;
    v_sc0 = orb_init.v_sc0;

    % Compute DSM true anomaly from the scanned angle
    offset = 1e-3;
    if apsis_flag == 1
        nu_m = pi + offset - nu_DSM;
    else
        nu_m = 2*pi + offset - nu_DSM;
    end

    % Wrap to [0, 2*pi)
    if nu_m < 0
        nu_m = nu_m + 2*pi;
    end

    % Time of flight from departure to DSM
    Em = 2.0 * atan(sqrt((1.0 - orb_init.ecc) / (1.0 + orb_init.ecc)) * tan(nu_m / 2.0));
    Mm = Em - orb_init.ecc * sin(Em);

    delta_M = Mm - orb_init.M0;
    if delta_M <= 0
        delta_M = delta_M + 2*pi;
    end
    delta_M = delta_M + (revs_before * 2*pi);
    
    tof_m_s = delta_M / orb_init.n_motion;

    % DSM must occur before arrival epoch
    if tof_m_s >= planets_state.tof_total_s
        return;
    end

    % Spacecraft state at DSM point
    [r_m_out, v_m_minus_out] = KEP2ICF_O(orb_init.sma, orb_init.ecc, orb_init.inc, nu_m, orb_init.argp, orb_init.raan, mu_sun);

    % Reject DSM points too close to the Sun
    AU = 149597870.7; % km
    if norm(r_m_out) < 0.2 * AU
        return;
    end

    % Remaining time of flight for Lambert arc
    tof_lambert_out = planets_state.tof_total_s - tof_m_s;

    % Multi-rev Lambert validity check
    mr_lambert = N - revs_before - 1;
    if mr_lambert > 0
        current_valid = checkLambertValidity(r_m_out, planets_state.r_pf, tof_lambert_out, mu_sun, mr_lambert);
    else
        current_valid = true;
    end
end
