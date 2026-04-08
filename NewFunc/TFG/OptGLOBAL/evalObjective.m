function [dV_total, dV_GA1, dV_DSM, dV_GA2, vinf_in, va, r_m, r_arc1, r_arc2, rp1_out, rp2_out] = evalObjective(vinf_out, vinfi, vinfin, vinfo, body, jd0, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, full_output_flag)
%   Unified objective function for evaluating the total Delta-V cost of a
%   V-Infinity Leveraging Maneuver (VILM) sequence dV_GA1 + dV_DSM + dV_GA2.
%     - Fast mode (full_output_flag = false): computes dV_total only, skipping arc generation.
%     - Full mode (full_output_flag = true): computes all outputs.
%
% Inputs:
%   vinf_out: outgoing v-infinity vector at the gravity assist body [km/s]
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
%   res_flag: resonance flag controlling GA1/GA2 cost
%   full_output_flag: boolean flag selecting fast or full output mode
%
% Outputs:
%   dV_total: total Delta-V cost of the VILM sequence [km/s]
%   dV_GA1: Delta-V cost at the first gravity assist (GA1) [km/s]
%   dV_DSM: Delta-V cost at the Deep Space Maneuver [km/s]
%   dV_GA2: Delta-V cost at the second gravity assist (GA2) [km/s]
%   vinf_in: inbound v-infinity vector at the resonant orbit arrival [km/s]
%   va: heliocentric arrival velocity of the spacecraft [km/s]
%   r_m: heliocentric position vector of the DSM [km]
%   r_arc1: position vectors along the pre-DSM trajectory arc [km]
%   r_arc2: position vectors along the post-DSM trajectory arc [km]
%   rp1_out: periapsis radius achieved at GA1 [km]
%   rp2_out: periapsis radius achieved at GA2 [km]
%
% References:
%   [-] n/a
%
% See also:
%   ResonantVILM2, OptimitzationVILM, GA_PGA2_Rp, GA_PGA2_Vinfo
%
% Adrià Solà Foixench
% April 2026

    % Default initialization of all outputs
    dV_total = inf; dV_GA1 = NaN; dV_DSM = NaN; dV_GA2 = NaN;
    vinf_in = [NaN NaN NaN]; va = [NaN NaN NaN]; r_m = [NaN NaN NaN];
    r_arc1 = []; r_arc2 = []; rp1_out = NaN; rp2_out = NaN;

    % Compute magnitude of the candidate outgoing v-infinity vector
    vinfon = norm(vinf_out);
    if vinfon < 1e-9
        return;
    end

    if full_output_flag
        % Full mode DSM Delta-V, arrival v-infinity, arrival velocity and trajectory arcs
        [dV_DSM, vinf_in, ~, va, r_m, r_arc1, r_arc2] = ResonantVILM2(body, jd0, vinf_out, N, M, apsis_flag, mu_sun, full_output_flag);
    else
        % Fast mode discard arc outputs
        [dV_DSM, vinf_in, ~, ~, ~, ~, ~] = ResonantVILM2(body, jd0, vinf_out, N, M, apsis_flag, mu_sun, full_output_flag);
    end

    if isnan(dV_DSM)
        return; 
    end

    % First Gravity Assist Delta-V cost
    if res_flag == 1
        dV_GA1 = 0;
    else
        % Compute deflection angle between vinfi and vinf_out
        if vinfin < 1e-9 
            return; 
        end
        cos_d1 = dot(vinfi, vinf_out) / (vinfin * vinfon);
        
        % Compute GA1 Delta-V and the corresponding periapsis radius
        [dV_GA1, rp1] = GA_PGA2_Rp(vinfin, vinfon, acos(cos_d1), mu_planet);
        if rp1 < vmr_safety
            % Periapsis below safety limit
            dV_GA1 = GA_PGA2_Vinfo(vinfin, vinfon, acos(cos_d1), vmr_safety, mu_planet);
            rp1_out = vmr_safety;
        else
            rp1_out = rp1;
        end
    end

    % Second Gravity Assist Delta-V cost
    if res_flag == 3
        dV_GA2 = 0;
    else
        vinf_in = vinf_in(:)';
        vinf_in_resn = norm(vinf_in);
        vinfon_req = norm(vinfo);
        
        if vinf_in_resn < 1e-9 || vinfon_req < 1e-9
            return; 
        end
        
        cos_d2 = dot(vinf_in, vinfo) / (vinf_in_resn * vinfon_req);
        
        % Compute GA2 Delta-V and the corresponding periapsis radius
        [dV_GA2, rp2] = GA_PGA2_Rp(vinf_in_resn, vinfon_req, acos(cos_d2), mu_planet);
        if rp2 < vmr_safety
            % Periapsis below safety limit
            dV_GA2 = GA_PGA2_Vinfo(vinf_in_resn, vinfon_req, acos(cos_d2), vmr_safety, mu_planet);
            rp2_out = vmr_safety;
        else
            rp2_out = rp2;
        end
    end

    % Accumulate total Delta-V cost across all three maneuver components
    dV_total = dV_GA1 + dV_DSM + dV_GA2;
end