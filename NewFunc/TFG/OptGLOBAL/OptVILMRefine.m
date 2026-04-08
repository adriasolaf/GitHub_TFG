function vinf_out_opt = OptVILMRefine(best_vinf_out, vinfi, vinfin, vinfo, body, jd0, N, M, ...
                                                apsis_flag, mu_sun, mu_planet, vmr_safety, vmag_min, vmag_max, res_flag)
%   Refinement of the VILM optimization. Takes the best outgoing v-infinity vector found during the phase 1 scan
%   and applies fminsearch to locate the nearest local minimum with higher precision.
%
% Inputs:
%   best_vinf_out: best Cartesian outgoing v-infinity vector from the initial scan [km/s]
%   vinfi: inbound v-infinity vector at the gravity assist body [km/s]
%   vinfin: magnitude of the inbound v-infinity vector [km/s]
%   vinfo: desired outbound v-infinity vector at GA2 [km/s]
%   body: identifier for the planetary body
%   jd0: initial departure epoch
%   N: number of spacecraft revolutions
%   M: number of planetary body revolutions
%   apsis_flag: flag defining apsis constraints for the DSM scan
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   mu_planet: planetary body gravitational parameter [km^3/s^2]
%   vmr_safety: minimum safe flyby radius with safety margin applied [km]
%   vmag_min: lower bound on the v-infinity magnitude search range [km/s]
%   vmag_max: upper bound on the v-infinity magnitude search range [km/s]
%   res_flag: resonance flag controlling GA1/GA2 cost computation
%
% Outputs:
%   vinf_out_opt: optimized outgoing v-infinity vector in Cartesian coordinates [km/s]
%
% References:
%   [-] n/a
%
% See also:
%   evalObjectiveFmin, evalObjective, OptVILMScan, OptVILMCompareAndGen, OptimitzationVILM
%
% Adrià Solà Foixench
% April 2026

    % Convert the best Cartesian v-infinity vector to spherical coordinates to form the initial guess for the optimizer
    [theta0, phi0, vmag0] = cart2sph(best_vinf_out(1), best_vinf_out(2), best_vinf_out(3));
    x0 = [vmag0, theta0, phi0];
    
    % Configure fminsearch convergence tolerances and iteration limit
    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 500, 'Display', 'off');
    
    obj_fun = @(x) evalObjectiveFmin(x, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, vmag_min, vmag_max, res_flag);
                                       
    x_opt = fminsearch(obj_fun, x0, options);
    
    % Reconstruct the optimal Cartesian v-infinity vector
    vmag_opt = max(vmag_min, min(vmag_max, x_opt(1)));
    theta_opt = x_opt(2);
    phi_opt = x_opt(3);
    
    vinf_out_opt = vmag_opt * [cos(phi_opt)*cos(theta_opt), cos(phi_opt)*sin(theta_opt), sin(phi_opt)];
end