function dV_total = objective_dV_total_fmin(x, vinfi, vinfo, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, vmag_min, vmag_max, res_flag, search_nu_handle)
%   Function that converts spherical coordinates (vmag, theta, phi) into
%   cartesian v-infinity vector and evaluates the total VILM cost with
%   objective_dV_total. Applies a penalty if the magnitude falls outside [vmag_min, vmag_max].
%
% Inputs:
%   x: optimization variable [vmag, theta, phi]
%          x(1) = v-infinity magnitude [km/s]
%          x(2) = azimuth angle theta [rad]
%          x(3) = elevation angle phi [rad]
%   vinfi: inbound v-infinity vector at GA1 [km/s]
%   vinfo: desired outbound v-infinity vector at GA2 [km/s]
%   body: identifier for the planetary body
%   jd0: departure epoch
%   T_p: orbital period of the resonant body [s]
%   N: number of spacecraft revolutions in the resonant orbit
%   M: number of resonant body revolutions in the resonant orbit
%   apsis_flag: DSM reference apsis (1 = apoapsis, 0 = periapsis)
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   mu_planet: planetary body gravitational parameter [km^3/s^2]
%   vmr_safety: minimum safe flyby radius [km]
%   vmag_min & vmag_max: bounds on v-infinity magnitude [km/s]
%   res_flag: resonance flag (1, 2, or 3)
%   search_nu_handle: handle to the 1D nu search function
%
% Outputs:
%   dV_total:  total Delta-V cost [km/s] (scalar, for fminsearch)
%
% References:
%   [-] n/a
%
% See also:
%   objective_dV_total, OptimitzationVILM
%
% Adria Sola Foixench
% April 2026

    vmag = x(1);
    theta = x(2);
    phi = x(3);

    % Convert spherical to Cartesian v-infinity vector
    vinf_out = vmag * [cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi)];

    % Evaluate total VILM cost in fast mode
    dV_total = objective_dV_total(vinf_out, vinfi, vinfo, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, search_nu_handle, false);

    % Penalize if magnitude is outside allowed bounds
    if vmag < vmag_min || vmag > vmag_max
        dV_total = dV_total + 1e6;
    end

end
