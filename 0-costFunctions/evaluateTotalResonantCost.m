function [dV_total, dV_GA1, dV_DSM, dV_GA2, vinf_in, va, r_m, v_m_minus, v_m_plus, rp1_out, rp2_out] = evaluateTotalResonantCost(vinf_out, vinfi_req, vinfo_req, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, search_nu_handle, full_output_flag)
%   Evaluates the total VILM Delta-V cost (dV_GA1 + dV_DSM + dV_GA2) for a
%   given outgoing v-infinity vector. Calls findOptimalDSMParameters to
%   find the optimal DSM, then computes the gravity assist costs.
%
% Inputs:
%   vinf_out: outgoing v-infinity vector at the resonant body [km/s]
%   vinfi_req: inbound v-infinity vector at GA1 [km/s]
%   vinfo_req: desired outbound v-infinity vector at GA2 [km/s]
%   body: identifier for the planetary body
%   jd0: departure epoch
%   T_p: orbital period of the resonant body [s]
%   N: number of spacecraft revolutions
%   M: number of body revolutions
%   apsis_flag: DSM reference apsis (1 = apoapsis, 0 = periapsis)
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   mu_planet: planetary body gravitational parameter [km^3/s^2]
%   vmr_safety: minimum safe flyby radius with safety margin [km]
%   res_flag: resonance flag:
%               1 = launch leg (dV_GA1 = 0)
%               2 = full VILM (both GAs computed)
%               3 = arrival leg (dV_GA2 = 0)
%   search_nu_handle: handle to the 1D nu search function
%   full_output_flag: true = return all vectors, false = return only costs
%
% Outputs:
%   dV_total: total Delta-V cost [km/s]
%   dV_GA1: Delta-V at the first gravity assist [km/s]
%   dV_DSM: Delta-V at the Deep Space Maneuver [km/s]
%   dV_GA2: Delta-V at the second gravity assist [km/s]
%   vinf_in: inbound v-infinity at the resonant orbit arrival [km/s]
%   va: heliocentric arrival velocity [km/s]
%   r_m: DSM position vector [km]
%   v_m_minus: velocity before DSM [km/s]
%   v_m_plus: velocity after DSM [km/s]
%   rp1_out: periapsis radius at GA1 [km]
%   rp2_out: periapsis radius at GA2 [km]
%
% References:
%   [-] n/a
%
% See also:
%   findOptimalDSMParameters, computeSingleDSMTransferCost, GA_PGA2_Rp, GA_PGA2_Vinfo,
%   optimizeOutgoingVInfinityVILM
%
% Adria Sola Foixench
% April 2026

    % Track how many times the cost is evaluated
    costEvalCounter('inc');

    % Initialize outputs
    dV_total = Inf;
    dV_GA1 = NaN;
    dV_DSM = NaN;
    dV_GA2 = NaN;
    vinf_in = [NaN NaN NaN];
    va = [NaN NaN NaN];
    r_m = [NaN NaN NaN];
    v_m_minus = [NaN NaN NaN];
    v_m_plus = [NaN NaN NaN];
    rp1_out = NaN;
    rp2_out = NaN;

    % Magnitudes of the v-infinity vectors
    vinfon = norm(vinf_out);
    vinfin = norm(vinfi_req);

    if vinfon < 1e-9
        return;
    end

    % 1. DSM cost
    [dV_DSM, ~, ~, ~, r_m_sol, v_m_minus_sol, v_m_plus_sol, vinf_in_sol, va_sol] = findOptimalDSMParameters(vinf_out, body, jd0, T_p, N, M, apsis_flag, mu_sun, search_nu_handle, true);

    if isnan(dV_DSM) || isinf(dV_DSM)
        dV_DSM = NaN;
        return;
    end

    
    % 2. GA1 cost: deflection from vinfi_req to vinf_out
    if res_flag == 1
        % Initial leg
        dV_GA1 = 0;
    else
        if vinfin < 1e-9
            return;
        end
        % Deflection angle between inbound and outgoing v-infinity

        cos_d1 = dot(vinfi_req, vinf_out) / (vinfin * vinfon);
        cos_d1 = max(-1.0, min(1.0, cos_d1));
        delta1 = acos(cos_d1);

        [dV_GA1, rp1] = GA_PGA2_Rp(vinfin, vinfon, delta1, mu_planet);

        % Enforce minimum flyby radius
        if rp1 < vmr_safety
            dV_GA1 = GA_PGA2_Vinfo(vinfin, vinfon, delta1, vmr_safety, mu_planet);
            rp1_out = vmr_safety;
        else
            rp1_out = rp1;
        end
    end

    
    % 3. GA2 cost: deflection from vinf_in to vinfo_req
    if res_flag == 3
        % Arrival leg
        dV_GA2 = 0;
    else
        vinf_in_vec = vinf_in_sol(:)';
        vinf_in_mag = norm(vinf_in_vec);
        vinfon_req = norm(vinfo_req);

        if vinf_in_mag < 1e-9 || vinfon_req < 1e-9
            return;
        end

        % Deflection angle between arrival v-infinity and required outbound
        cos_d2 = dot(vinf_in_vec, vinfo_req) / (vinf_in_mag * vinfon_req);
        cos_d2 = max(-1.0, min(1.0, cos_d2));
        delta2 = acos(cos_d2);

        [dV_GA2, rp2] = GA_PGA2_Rp(vinf_in_mag, vinfon_req, delta2, mu_planet);

        % Enforce minimum flyby radius
        if rp2 < vmr_safety
            dV_GA2 = GA_PGA2_Vinfo(vinf_in_mag, vinfon_req, delta2, vmr_safety, mu_planet);
            rp2_out = vmr_safety;
        else
            rp2_out = rp2;
        end
    end

    % 4. Total cost
    dV_total = dV_GA1 + dV_DSM + dV_GA2;

    % Fill output vectors if requested
    if full_output_flag
        vinf_in = vinf_in_sol;
        va = va_sol;
        r_m = r_m_sol;
        v_m_minus = v_m_minus_sol;
        v_m_plus = v_m_plus_sol;
    end

end
