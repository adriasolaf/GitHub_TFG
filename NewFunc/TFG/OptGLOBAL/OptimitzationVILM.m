function [best_vinf_out, dV_total_nave, best_dV_GA1, best_dV_DSM, best_dV_GA2, best_vinf_in, best_va, best_vinf_o_n, orbita_res, rp_GA1, rp_GA2] = OptimitzationVILM(body, jd0, vinfi, vinfo, N, M, apsis_flag, mu_sun, res_flag)
%   Optimizes the outgoing v-infinity vector at a planetary flyby to minimize
%   the total Delta-V cost of an N:M resonant maneuver.
%   It optimizes the full VILM cost: dV_GA1 + dV_DSM + dV_GA2.
%   The DSM Delta-V minimization for a given vinf_out is done by ResonantVILM2.
%
% Inputs:
%   body: identifier for the planetary body
%   jd0: initial departure epoch
%   vinfi: inbound v-infinity vector at the gravity assist body [km/s]
%   vinfo: desired outbound v-infinity vector at the gravity assist body [km/s]
%   N: number of spacecraft revolutions
%   M: number of planetary body revolutions
%   apsis_flag: flag defining apsis constraints for the DSM scan
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   res_flag: resonance flag controlling GA1/GA2 cost:
%                       1 = launch leg (GA1 cost = 0)
%                       3 = arrival leg (GA2 cost = 0)
%                       otherwise = full two-flyby VILM sequence
%
% Outputs:
%   best_vinf_out: optimal outgoing v-infinity vector at departure [km/s]
%   dV_total_nave: total Delta-V cost of the optimal VILM sequence [km/s]
%   best_dV_GA1: Delta-V cost at the first gravity assist [km/s]
%   best_dV_DSM: Delta-V cost at the Deep Space Maneuver [km/s]
%   best_dV_GA2: Delta-V cost at the second gravity assist (GA2) [km/s]
%   best_vinf_in: optimal v-infinity vector at arrival [km/s]
%   best_va: optimal heliocentric arrival velocity of the spacecraft [km/s]
%   best_vinf_o_n: magnitude of the optimal outgoing v-infinity at GA2 [km/s]
%   orbita_res: resonant orbit parameters struct
%   rp_GA1: periapsis radius at the first gravity assist [km]
%   rp_GA2: periapsis radius at the second gravity assist [km]
%
% References:
%   [-] n/a
%
% See also:
%   ResonantVILM2, GetBodyProps, OptVILMScan, OptVILMRefine, OptVILMCompareAndGen
%
% Adrià Solà Foixench
% April 2026


    % Planetary gravitational parameter and minimum flyby radius
    [mu_planet, vmr] = GetBodyProps(body);
    % Apply safety margin to minimum flyby radius to avoid atmospheric entry
    vmr_safety = 1.05 * vmr;
    % Compute scalar magnitudes of inbound and outbound v-infinity vectors
    vinfin = norm(vinfi);
    vinfon_required = norm(vinfo);

    if res_flag == 1
        % Only outbound v-infinity is defined
        vmag_min = max(0.1, vinfon_required * 0.3);
        vmag_max = max(5.0, vinfon_required * 1.5);
        v_ref = vinfo; % Uses vinfo as reference since vinfi doesn't exist
    elseif res_flag == 3
        % Only inbound v-infinity is defined
        vmag_min = max(0.1, vinfin * 0.3);
        vmag_max = max(5.0, vinfin * 1.5);
        v_ref = vinfi; % Uses vinfi as reference
    else
        % Both v-infinity vectors are defined
        vmag_min = max(0.1, min(vinfin, vinfon_required) * 0.3);
        vmag_max = max(vinfin, vinfon_required) * 1.5;
        v_ref = vinfi;
    end

    % Number of magnitude, azimuth and elevation samples
    N_mag = 10;
    N_theta = 24;
    N_phi = 5;

    % Normalize the reference direction vector
    v_ref_n = norm(v_ref);
    if v_ref_n > 1e-9
        v_ref_dir = v_ref / v_ref_n;
    else
        v_ref_dir = [1, 0, 0]; % Default for nul cases
    end
    % Convert reference direction to spherical coordinates
    [theta_ref, phi_ref, ~] = cart2sph(v_ref_dir(1), v_ref_dir(2), v_ref_dir(3));

    % Check if the reference vector is in the ecliptic plane
    in_plane = (abs(v_ref_dir(3))) < 1e-6;

    % Build magnitude and azimuth search grids centered on the reference direction
    vmag_vec = linspace(vmag_min, vmag_max, N_mag);
    theta_vec = linspace(theta_ref - pi, theta_ref + pi, N_theta + 1);
    theta_vec = theta_vec(1:end-1); % Remove duplicate

    if in_plane
        phi_vec = 0;
    else
        phi_vec = linspace(phi_ref - pi/4, phi_ref + pi/4, N_phi);
    end

    % Total number of magnitude samples
    total_vmag = length(vmag_vec);
    
    % Evaluate VILM cost over vmag, theta and phi
    [best_dV_total_scan, best_vinf_out_scan] = OptVILMScan(vmag_vec, total_vmag, theta_vec, phi_vec, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag);

    if isinf(best_dV_total_scan)
        return;
    end

    % Optimize vinf_out starting from the best solution
    vinf_out_opt = OptVILMRefine(best_vinf_out_scan, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, vmag_min, vmag_max, res_flag);

    % Select the global best and generate outputs
    [best_vinf_out, dV_total_nave, best_dV_GA1, best_dV_DSM, best_dV_GA2, best_vinf_in, best_va, best_vinf_o_n, orbita_res, rp_GA1, rp_GA2] = OptVILMCompareAndGen(vinf_out_opt, best_vinf_out_scan, best_dV_total_scan, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag);
end