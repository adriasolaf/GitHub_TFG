function [dV, vinf_in, v_sc0, va, r_m_out, r_arc1, r_arc2] = ResonantVILM2(body, jd0, vinf_out, N, M, apsis_flag, mu_sun, full_output_flag)
%   Searches for the optimal DSM location and timing for an N:M resonant maneuver. 
%   The optimization criterion is the minimum of the DSM Delta-V. The optimization of 
%   dV_GA1 + dV_DSM + dV_GA2 is done in OptimitzationVILM.
%
% Inputs:
%   body: identifier for the planetary body
%   jd0: initial departure epoch
%   vinf_out: v-infinity vector at departure [km/s]
%   N: number of spacecraft revolutions
%   M: number of planetary body revolutions
%   apsis_flag: flag defining apsis constraints for the scan
%   mu_sun: central body gravitational parameter [km^3/s^2]
%
% Outputs:
%   dV: magnitude of the optimal DSM Delta-V [km/s]
%   vinf_in: inbound v-infinity vector at arrival [km/s]
%   v_sc0: initial heliocentric velocity of the s/c at departure [km/s]
%   va: heliocentric velocity of the s/c at arrival [km/s]
%   r_m_out: heliocentric position vector of the optimal DSM [km]
%   r_arc1: position vectors along the pre-DSM trajectory arc [km]
%   r_arc2: position vectors along the post-DSM trajectory arc [km]
%
% References:
%   [1] PhD_DavidDLTS_astro (David de la Torre) [2] Trajectory design of solar orbiter (José
%   Manuel Sánchez Pérez)
%
% See also:
%   OptimitzationVILM, GetBodyICF, ICF2KEP_O, scan_nu_pos, find_best_lambert
%
% Adrià Solà Foixench
% April 2026
 
% Constants
sec2days = 1 / 86400; % Seconds to days conversion factor
T_p = 11.86 * 365.25636 * 24 * 3600; % Jupiter orbital period [s]

% Planet position and velocity at departure epoch
[r_p0, v_p0] = GetBodyICF(body, jd0, mu_sun, 0);
 
% Total time of flight
tof_total_s = M * T_p;
jd_f = jd0 + (tof_total_s * sec2days); % Final epoch
 
% Planet position and velocity at arrival epoch
[r_pf, v_pf] = GetBodyICF(body, jd_f, mu_sun, 0);

% Spacecraft initial orbit
v_sc0 = v_p0 + vinf_out; % Spacecraft heliocentric velocity at departure [km/s]
 
% Keplerian elements of the initial heliocentric orbit
[sma, ecc, inc, nu0, argp, raan] = ICF2KEP_O(r_p0, v_sc0, mu_sun);
 
% Check for hyperbolic orbit
if ecc >= 1.0
    dV = NaN; vinf_in = [NaN NaN NaN]; va = [NaN NaN NaN];
    r_m_out = [NaN NaN NaN]; r_arc1 = []; r_arc2 = [];
    return;
end
 
% Mean motion and initial mean anomaly
n = sqrt(mu_sun / (sma^3)); % Mean motion [rad/s]
E0 = 2.0 * atan(sqrt((1.0 - ecc) / (1.0 + ecc)) * tan(nu0 / 2.0)); % Eccentric anomaly [rad]
M0 = E0 - ecc * sin(E0); % Mean anomaly [rad]
 
% Change angle limits for DSM search
offset = 1e-3; % Small offset to avoid singularities at apsis
 
% Maximum change angle
change_angle_max = pi - offset; % Search up to pi

% Search grid parameters
N_points = 500; % Number of points
N_refine = 60; % Number of refinement iterations

% Initialize best solution variables
best_dV_global = inf; % Best dV
revs_before_m = 0; % Revolutions before DSM point
mr_lambert = 0; % Multi-revolution Lambert count
r_m = zeros(1,3); % DSM position [km]
v_m_minus = zeros(1,3); % Velocity before DSM [km/s]
dV = NaN; % Optimal DSM Delta-V magnitude [km/s]
vinf_in = [NaN NaN NaN]; % Inbound v-infinity vector at arrival [km/s]
va = [NaN NaN NaN]; % Heliocentric arrival velocity of the spacecraft [km/s]
v_m_plus = [NaN NaN NaN]; % Spacecraft velocity immediately after the DSM [km/s]

% Initialize structs
orb.sma = sma; orb.ecc = ecc; orb.inc = inc; orb.nu0 = nu0; orb.argp = argp; orb.raan = raan; orb.n = n; orb.M0 = M0;
aux.apsis_flag = apsis_flag; aux.offset = offset; aux.r_pf = r_pf; aux.tof_total_s = tof_total_s; aux.jd0 = jd0; aux.jd_f = jd_f; aux.sec2days = sec2days; aux.mu_sun = mu_sun; 
 
% Loop over possible revolution splits
for revs_try = N-1:-1:0
 
    mr_lambert_try = N - revs_try - 1; % Lambert multi-rev
    lp_max_try = 0; % Long-period flag
    if mr_lambert_try > 0
        lp_max_try = 1; % Long-period of the two posible branches
    end
 
    nu_opt = scan_nu_pos(N_points, N_refine, change_angle_max, orb, aux, revs_try, mr_lambert_try, lp_max_try);

    if isnan(nu_opt)
        continue; % Skip if no valid solution found in this revs_try
    end

    % Evaluate final solution at refined point
    [valid_opt, r_m_opt, v_m_minus_opt, tof_opt] = checkValidNu(nu_opt, orb, aux, revs_try, mr_lambert_try);
 
    if ~valid_opt
        continue; 
    end
   
    % Try all Lambert solution branches
    [dV_opt, vinf_in_opt, va_opt, v_m_plus_opt] = find_best_lambert(r_m_opt, r_pf, tof_opt, mu_sun, v_m_minus_opt, v_pf, mr_lambert_try, lp_max_try);
 
    % Skip if no valid Lambert
    if isinf(dV_opt)
        continue;
    end
 
    % Update global best solution
    if dV_opt < best_dV_global
        best_dV_global = dV_opt;
        revs_before_m = revs_try;
        mr_lambert = mr_lambert_try;
        r_m = r_m_opt;
        v_m_minus = v_m_minus_opt;
        dV = dV_opt;
        vinf_in = vinf_in_opt;
        va = va_opt;
        v_m_plus = v_m_plus_opt;
    end
 
end 
 
% Check if any valid solution was found
if isinf(best_dV_global)
    dV = NaN; vinf_in = [NaN NaN NaN]; va = [NaN NaN NaN];
    r_m_out = [NaN NaN NaN]; r_arc1 = []; r_arc2 = [];
    return;
end

% Compute position of resonant arcs
if full_output_flag
    [r_arc1, r_arc2, r_m_out] = Resonant_arcs(r_p0, v_sc0, r_m, v_m_minus, v_m_plus, r_pf, va, revs_before_m, mr_lambert, mu_sun);
else 
    r_arc1 = []; r_arc2 = []; r_m_out = r_m;
end
end 
 
 
