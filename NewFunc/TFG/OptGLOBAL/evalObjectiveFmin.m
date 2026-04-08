function dV_total = evalObjectiveFmin(x, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, vmag_min, vmag_max, res_flag)
%   Converts the spherical vmag, theta, phi into a Cartesian
%   outgoing v-infinity vector and evaluates the total VILM Delta-V cost via evalObjective.
%
% Inputs:
%   x: optimization variable vector [vmag, theta, phi]
%                  x(1): v-infinity magnitude [km/s]
%                  x(2): azimuth angle theta [rad]
%                  x(3): elevation angle phi [rad]
%   vinfi: inbound v-infinity vector at the gravity assist body [km/s]
%   vinfin: magnitude of the inbound v-infinity vector [km/s]
%   vinfo: desired outbound v-infinity vector at GA2 [km/s]
%   body: identifier for the planetary body
%   jd0: initial departure epoch [Julian Date]
%   N: number of spacecraft revolutions
%   M: number of planetary body revolutions
%   apsis_flag: flag defining apsis constraints for the DSM scan
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   mu_planet: planetary body gravitational parameter [km^3/s^2]
%   vmr_safety: minimum safe flyby radius with safety margin applied [km]
%   vmag_min: lower bound on the v-infinity magnitude search range [km/s]
%   vmag_max: upper bound on the v-infinity magnitude search range [km/s]
%   res_flag: resonance flag GA1/GA2 cost computation
%
% Outputs:
%   dV_total: total Delta-V cost of the VILM sequence [km/s]
%
% References:
%   [-] n/a
%
% See also:
%   evalObjective, OptVILMRefine, OptimitzationVILM
%
% Adrià Solà Foixench
% April 2026

    vmag = x(1);
    theta = x(2);
    phi = x(3);

    vinf_out = vmag * [cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi)];

    % Evaluate the total VILM Delta-V cost
    dV_total = evalObjective(vinf_out, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, false);

    if x(1) < vmag_min || x(1) > vmag_max
        dV_total = dV_total + 1e6; % Penalty
    end
end