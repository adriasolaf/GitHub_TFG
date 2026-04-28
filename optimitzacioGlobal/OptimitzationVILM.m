function [best_vinf_out, dV_total, best_dV_GA1, best_dV_DSM, best_dV_GA2, best_vinf_in, best_va, best_vinf_o_n, orbit_res, rp_GA1, rp_GA2] = OptimitzationVILM(body, jd0, vinfi, vinfo, N, M, apsis_flag, mu_sun, res_flag)
%   Optimizes the outgoing v-infinity vector at a planetary flyby to minimize
%   the total Delta-V cost of an N:M resonant VILM maneuver: dV_total = dV_GA1 + dV_DSM + dV_GA2
%   Strategy: 3D grid scan in spherical coordinates (vmag, theta, phi)
%             to locate the global minimum rapidly. Then, local refinement with fminsearch starting from the
%             previous best point. Finally, comparison of refined vs scan solution.
%
% Inputs:
%   body: identifier for the planetary body
%   jd0: departure epoch
%   vinfi: inbound v-infinity vector at GA1 [km/s]
%   vinfo: desired outbound v-infinity vector at GA2 [km/s]
%   N: number of spacecraft revolutions
%   M: number of body revolutions
%   apsis_flag: DSM reference apsis (1 = apoapsis, 0 = periapsis)
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   res_flag: resonance flag (1 = departure, 2 = intermediate, 3 = arrival)
%
% Outputs:
%   best_vinf_out: optimal outgoing v-infinity vector [km/s]
%   dV_total: total Delta-V cost [km/s]
%   best_dV_GA1: Delta-V at GA1 [km/s]
%   best_dV_DSM: Delta-V at DSM [km/s]
%   best_dV_GA2: Delta-V at GA2 [km/s]
%   best_vinf_in: arrival v-infinity vector [km/s]
%   best_va: arrival heliocentric velocity [km/s]
%   best_vinf_o_n: magnitude of the optimal outgoing v-infinity [km/s]
%   orbit_res: struct with trajectory arc data (.rdsm, .bdsm, .adsm)
%   rp_GA1: periapsis radius at GA1 [km]
%   rp_GA2: periapsis radius at GA2 [km]
%
% References:
%   [-] n/a
%
% See also:
%   objective_dV_total, objective_dV_total_fmin, solve_best_DSM,
%   grid_scan_3d, scan_and_refine_1d
%
% Adria Sola Foixench
% April 2026

    % Planetary parameters
    [mu_planet, vmr] = GetBodyProps(body);
    vmr_safety = 1.05 * vmr; % Safety margin on minimum flyby radius

    % Orbital period of the resonant body
    [sma_planet, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(body, jd0);
    T_p = 2*pi*sqrt(sma_planet^3 / mu_sun);

   
    % Search bounds for vinf_out based on available v-infinity magnitudes
    vinfin = norm(vinfi);
    vinfon_required = norm(vinfo);

    if res_flag == 1
        vmag_min = max(0.1, vinfon_required * 0.3);
        vmag_max = max(5.0, vinfon_required * 1.5);
        v_ref = vinfo;
    elseif res_flag == 3
        vmag_min = max(0.1, vinfin * 0.3);
        vmag_max = max(5.0, vinfin * 1.5);
        v_ref = vinfi;
    else
        vmag_min = max(0.1, min(vinfin, vinfon_required) * 0.3);
        vmag_max = max(vinfin, vinfon_required) * 1.5;
        v_ref = vinfi;
    end

   
    % Reference direction in spherical coordinates
    v_ref_n = norm(v_ref);
    if v_ref_n > 1e-9
        v_ref_dir = v_ref / v_ref_n;
    else
        v_ref_dir = [1, 0, 0];
    end
    [theta_ref, phi_ref, ~] = cart2sph(v_ref_dir(1), v_ref_dir(2), v_ref_dir(3));

    % If the reference vector is in the ecliptic plane, skip phi search
    in_plane = (abs(v_ref_dir(3))) < 1e-6;

    
    % Build search grids
    N_mag = 10;
    N_theta = 24;
    N_phi = 5;

    vmag_vec = linspace(vmag_min, vmag_max, N_mag);
    theta_vec = linspace(theta_ref - pi, theta_ref + pi, N_theta + 1);
    theta_vec = theta_vec(1:end-1); % Remove duplicate endpoint

    if in_plane
        phi_vec = 0;
    else
        phi_vec = linspace(phi_ref - pi/4, phi_ref + pi/4, N_phi);
    end

    % Define the 1D nu search method
    N_points_nu = 500;
    N_refine_nu = 60;
    search_nu = @(f, a, b) scan_and_refine_1d(f, a, b, N_points_nu, N_refine_nu);

   
    % 3D grid scan
    % Build objective for the 3D scan: x = [vmag, theta, phi]
    f_scan = @(x) objective_dV_total(x(1) * [cos(x(3))*cos(x(2)), cos(x(3))*sin(x(2)), sin(x(3))], vinfi, vinfo, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, search_nu, false);

    [x_best_scan, dV_best_scan] = grid_scan_3d(f_scan, vmag_vec, theta_vec, phi_vec);

    % Reconstruct best cartesian vector from scan
    vinf_out_scan = x_best_scan(1) * [cos(x_best_scan(3))*cos(x_best_scan(2)), cos(x_best_scan(3))*sin(x_best_scan(2)), sin(x_best_scan(3))];

    if isinf(dV_best_scan)
        % No valid solution found in the entire grid
        best_vinf_out = [NaN NaN NaN];
        dV_total = Inf;
        best_dV_GA1 = NaN;
        best_dV_DSM = NaN;
        best_dV_GA2 = NaN;
        best_vinf_in = [NaN NaN NaN];
        best_va = [NaN NaN NaN];
        best_vinf_o_n = NaN;
        orbit_res = struct();
        rp_GA1 = NaN;
        rp_GA2 = NaN;
        return;
    end

    
    % Local refinement with fminsearch
    [theta0, phi0, vmag0] = cart2sph(vinf_out_scan(1), vinf_out_scan(2), vinf_out_scan(3));

    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 500, 'Display', 'off');

    if in_plane
        % 2D: optimize only [vmag, theta]
        obj_fmin = @(x) objective_dV_total_fmin([x(1), x(2), 0], vinfi, vinfo, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, vmag_min, vmag_max, res_flag, search_nu);
        x_opt2d = fminsearch(obj_fmin, [vmag0, theta0], options);
        x_opt = [x_opt2d(1), x_opt2d(2), 0];
    else
        obj_fmin = @(x) objective_dV_total_fmin(x, vinfi, vinfo, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, vmag_min, vmag_max, res_flag, search_nu);
        x0 = [vmag0, theta0, phi0];
        x_opt = fminsearch(obj_fmin, x0, options);
    end

    % Reconstruct refined cartesian vector
    vmag_opt = max(vmag_min, min(vmag_max, x_opt(1)));
    theta_opt = x_opt(2);
    phi_opt = x_opt(3);
    vinf_out_refined = vmag_opt * [cos(phi_opt)*cos(theta_opt), cos(phi_opt)*sin(theta_opt), sin(phi_opt)];

   
    % Phase 3: Compare refined vs scan and select best
    % Evaluate refined solution in full mode
    [dV_refined, dV_GA1_r, dV_DSM_r, dV_GA2_r, vinf_in_r, va_r, ~, ~, ~, rp1_r, rp2_r] = objective_dV_total(vinf_out_refined, vinfi, vinfo, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, search_nu, true);

    if dV_refined < dV_best_scan
        % Refined solution is better
        best_vinf_out = vinf_out_refined;
        dV_total = dV_refined;
        best_dV_GA1 = dV_GA1_r;
        best_dV_DSM = dV_DSM_r;
        best_dV_GA2 = dV_GA2_r;
        best_vinf_in = vinf_in_r;
        best_va = va_r;
        rp_GA1 = rp1_r;
        rp_GA2 = rp2_r;
    else
        % Scan solution is better; re-evaluate in full mode
        [dV_scan_full, dV_GA1_s, dV_DSM_s, dV_GA2_s, vinf_in_s, va_s, ~, ~, ~, rp1_s, rp2_s] = objective_dV_total(vinf_out_scan, vinfi, vinfo, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, search_nu, true);

        best_vinf_out = vinf_out_scan;
        dV_total = dV_scan_full;
        best_dV_GA1 = dV_GA1_s;
        best_dV_DSM = dV_DSM_s;
        best_dV_GA2 = dV_GA2_s;
        best_vinf_in = vinf_in_s;
        best_va = va_s;
        rp_GA1 = rp1_s;
        rp_GA2 = rp2_s;
    end

    % Build trajectory arcs for visualization
    best_vinf_o_n = norm(best_vinf_out);

    % Get revs_best and lp/lw from the solver to pass to Resonant_arcs
    [~, ~, revs_best, ~, ~, r_m_arc, v_m_minus_arc, v_m_plus_arc, ~, va_arc] = solve_best_DSM(best_vinf_out, body, jd0, T_p, N, M, apsis_flag, mu_sun, search_nu, true);

    mr_lambert_best = N - revs_best - 1;

    % Departure state
    [r_p0, v_p0] = GetBodyICF(body, jd0, mu_sun, 0);
    v_sc0 = v_p0 + best_vinf_out;

    % Arrival planet position
    sec2days = 1 / 86400;
    tof_total_s = M * T_p;
    jd_f = jd0 + (tof_total_s * sec2days);
    [r_pf, ~] = GetBodyICF(body, jd_f, mu_sun, 0);

    [r_arc1, r_arc2, rdsm] = Resonant_arcs(r_p0, v_sc0, r_m_arc, v_m_minus_arc, v_m_plus_arc, r_pf, va_arc, revs_best, mr_lambert_best, mu_sun);

    orbit_res.rdsm = rdsm;
    orbit_res.bdsm = r_arc1;
    orbit_res.adsm = r_arc2;

end
