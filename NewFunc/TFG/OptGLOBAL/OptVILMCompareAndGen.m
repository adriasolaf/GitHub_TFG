function [best_vinf_out, dV_total_nave, best_dV_GA1, best_dV_DSM, best_dV_GA2, best_vinf_in, best_va, best_vinf_o_n, orbita_res, rp_GA1, rp_GA2] = OptVILMCompareAndGen(vinf_out_opt, best_vinf_out_scan, best_dV_total_scan, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag)
%   Compares the refined solution against the best solution of phase 1, selects the global minimum,
%   and generates the full VILM output.
% Inputs:
%   vinf_out_opt: refined outgoing v-infinity vector from the local optimizer [km/s]
%   best_vinf_out_scan: best outgoing v-infinity vector from the phase 1 scan [km/s]
%   best_dV_total_scan: total Delta-V of the best phase 1 solution [km/s]
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
%   best_vinf_out: optimal outgoing v-infinity vector [km/s]
%   dV_total_nave: total Delta-V cost of the selected optimal VILM sequence [km/s]
%   best_dV_GA1: Delta-V cost at the first gravity assist (GA1) [km/s]
%   best_dV_DSM: Delta-V cost at the Deep Space Maneuver [km/s]
%   best_dV_GA2: Delta-V cost at the second gravity assist (GA2) [km/s]
%   best_vinf_in: inbound v-infinity vector at the resonant orbit arrival [km/s]
%   best_va: heliocentric arrival velocity of the spacecraft [km/s]
%   best_vinf_o_n: magnitude of the optimal outgoing v-infinity vector [km/s]
%   orbita_res: struct containing the resonant trajectory arc data:
%                       .rdsm - DSM heliocentric position vector [km]
%                       .bdsm - position vectors of the pre-DSM arc [km]
%                       .adsm - position vectors of the post-DSM arc [km]
%   rp_GA1: periapsis radius achieved at GA1 [km]
%   rp_GA2: periapsis radius achieved at GA2 [km]
%
% References:
%   [-] n/a
%
% See also:
%   evalObjective, OptVILMRefine, OptVILMScan, OptimitzationVILM
%
% Adrià Solà Foixench
% April 2026

    % Re-evaluate the refined solution in full mode to obtain all outputs including arcs
    [dV_tot_final, dV_GA1_final, dV_DSM_final, dV_GA2_final, vinf_in_final, va_final, r_m_final, arc1_final, arc2_final, rp1_final, rp2_final] = evalObjective(vinf_out_opt, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, true);

    % Accept the refined solution only if it improves result
    if dV_tot_final < best_dV_total_scan
        % Store all refined outputs
        best_dV_total = dV_tot_final;
        best_vinf_out = vinf_out_opt;
        best_dV_GA1 = dV_GA1_final;
        best_dV_DSM = dV_DSM_final;
        best_dV_GA2 = dV_GA2_final;
        best_vinf_in = vinf_in_final;
        best_va = va_final;
        rdsm_save = r_m_final;
        bdsm = arc1_final;
        adsm = arc2_final;
        rp_GA1 = rp1_final;
        rp_GA2 = rp2_final;
    else
        % Back to the phase 1 best point and evaluate in full mode
        [dV_check, best_dV_GA1, best_dV_DSM, best_dV_GA2, best_vinf_in, best_va, rdsm_save, bdsm, adsm, rp_GA1, rp_GA2] = evalObjective(best_vinf_out_scan, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, true);
                          
        if isinf(dV_check)
            best_dV_total = inf;
        else
            best_dV_total = dV_check;
        end
        best_vinf_out = best_vinf_out_scan;
    end

    % Assign final outputs and trajectory arc data into a struct
    dV_total_nave = best_dV_total;
    best_vinf_o_n = norm(best_vinf_out);
    orbita_res.rdsm = rdsm_save;
    orbita_res.bdsm = bdsm;
    orbita_res.adsm = adsm;
end