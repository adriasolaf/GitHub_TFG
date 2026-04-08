function [best_dV_total, best_vinf_out] = OptVILMScan(vmag_vec, total_vmag, theta_vec, phi_vec, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag)
%   Evaluates the total VILM Delta-V cost of outgoing v-infinity vectors, defined by magnitude (vmag),
%   theta and phi samples. The global minimum found is returned as the starting point for the refinement stage. 
%   All evaluations use the fast mode of evalObjective to minimize computation time.
%
% Inputs:
%   vmag_vec: vector of v-infinity magnitude samples [km/s]
%   total_vmag: number of magnitude samples (length of vmag_vec)
%   theta_vec: vector of azimuth angle samples [rad]
%   phi_vec: vector of elevation angle samples [rad]
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
%   res_flag: resonance flag controlling GA1/GA2 cost computation
%
% Outputs:
%   best_dV_total: minimum total Delta-V found [km/s]
%   best_vinf_out: Cartesian outgoing v-infinity vector achieving the minimum [km/s]
%
% References:
%   [-] n/a
%
% See also:
%   evalObjective, OptVILMRefine, OptVILMCompareAndGen, OptimitzationVILM
%
% Adrià Solà Foixench
% April 2026

    % Initialize global best values to infeasible defaults
    best_dV_total = inf;
    best_vinf_out = [NaN NaN NaN];
    
    % Sweep over all vmag, theta, phi combinations
    for idx_v = 1:total_vmag
        vmag = vmag_vec(idx_v);
        
        for theta = theta_vec
            for phi = phi_vec
                % Convert spherical grid point to a Cartesian candidate v-infinity vector
                vinf_out_cand = vmag * [cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi)];
                
                % Evaluate total VILM cost in fast mode
                [dV_total_cand, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = evalObjective(vinf_out_cand, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, false);
                                    
                % Update global best if this candidate improves the current minimum
                if dV_total_cand < best_dV_total
                    best_dV_total = dV_total_cand;
                    best_vinf_out = vinf_out_cand;
                end
            end
        end
    end
    
end