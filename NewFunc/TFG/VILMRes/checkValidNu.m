function [current_valid, r_m_out, v_m_minus_out, tof_lambert_out] = checkValidNu(change_angle, orb, aux, revs_try, mr_lambert_try)
%   Computes the true anomaly, calculates the time of flight to this 
%   DSM point, and determines if the resulting Lambert transfer to the arrival planet is valid.
%
% Inputs:
%   change_angle: the scanned angle used to determine the true anomaly [rad]
%   orb: structure containing the initial heliocentric Keplerian elements 
%        of the spacecraft
%   aux: structure containing auxiliary simulation parameters
%   revs_try: number of full spacecraft revolutions before the DSM point
%   mr_lambert_try: number of full revolutions for the Lambert transfer arc
%
% Outputs:
%   current_valid: boolean flag indicating if the DSM is valid
%   r_m_out: heliocentric position vector at the DSM point [km]
%   v_m_minus_out: heliocentric velocity vector just before the DSM [km/s]
%   tof_lambert_out: remaining time of flight for the Lambert arc [s]
%
% References:
%   [-] n/a
%
% See also:
%   checkLambertValidity, scan_nu_pos
%
% Adrià Solà Foixench
% April 2026

    % Initialize outputs
    current_valid = false;
    r_m_out = zeros(1,3);
    v_m_minus_out = zeros(1,3);
    tof_lambert_out = 0;

    % Compute true anomaly at DSM from nu angle
    if aux.apsis_flag == 1
        nu_m = pi + aux.offset - change_angle; % Periapsis reference
    else
        nu_m = 2*pi + aux.offset - change_angle; % Apoapsis reference
    end

    % nu [0, 2*pi)
    if nu_m < 0
        nu_m = nu_m + 2*pi;
    end

    % Compute time of flight to DSM point
    % Eccentric anomaly at DSM
    Em = 2.0 * atan(sqrt((1.0 - orb.ecc) / (1.0 + orb.ecc)) * tan(nu_m / 2.0));
    % Mean anomaly at DSM
    Mm = Em - orb.ecc * sin(Em);

    % Mean anomaly difference with revolution count
    delta_M = Mm - orb.M0;
    if delta_M <= 0
        delta_M = delta_M + 2*pi;
    end
    delta_M = delta_M + (revs_try * 2*pi);

    % Time of flight to DSM [s] and DSM epoch
    tof_m_s = delta_M / orb.n;
    jd_m = aux.jd0 + (tof_m_s * aux.sec2days);

    % DSM must occur before arrival
    if jd_m >= aux.jd_f
        return;
    end

    % Compute DSM position and velocity on the initial orbit
    [r_m_out, v_m_minus_out] = KEP2ICF_O(orb.sma, orb.ecc, orb.inc, nu_m, orb.argp, orb.raan, aux.mu_sun);

    % Remaining time of flight for Lambert arc [s]
    tof_lambert_out = aux.tof_total_s - tof_m_s;

    % Additional checks for multi-revolution Lambert
    if mr_lambert_try > 0
        current_valid = checkLambertValidity(r_m_out, aux.r_pf, tof_lambert_out, aux.mu_sun, mr_lambert_try);
    else
        current_valid = true; % Always valid if time checks pass
    end
end