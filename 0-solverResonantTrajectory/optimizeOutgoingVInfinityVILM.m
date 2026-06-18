function [best_vinf_out, dV_total, best_dV_GA1, best_dV_DSM, best_dV_GA2, best_vinf_in, best_va, best_vinf_o_n, vilm_arc, rp_GA1, rp_GA2] = optimizeOutgoingVInfinityVILM(body, jd0, vinfi_req, vinfo_req, N, M, apsis_flag, mu_sun, res_flag, opts)
%   Optimizes the outgoing v-infinity vector at a planetary flyby to minimize
%   the total Delta-V cost of an N:M resonant VILM maneuver: dV_total = dV_GA1 + dV_DSM + dV_GA2
%   Strategy: 3D grid scan in spherical coordinates (vmag, theta, phi)
%             to locate the global minimum rapidly. Then, optional local
%             refinement with fminsearch starting from the previous best
%             point, comparing refined vs scan solution.
%
% Inputs:
%   body: identifier for the planetary body
%   jd0: departure epoch
%   vinfi_req: inbound v-infinity vector at GA1 [km/s]
%   vinfo_req: desired outbound v-infinity vector at GA2 [km/s]
%   N: number of spacecraft revolutions
%   M: number of body revolutions
%   apsis_flag: DSM reference apsis (1 = apoapsis, 0 = periapsis)
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   res_flag: resonance flag (1 = departure, 2 = intermediate, 3 = arrival)
%   opts: struct with optimizer settings.
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
%   vilm_arc: struct packing the inputs needed to rebuild the resonant
%             trajectory arcs externally and to replay the optimizer in
%             validation tests (vinf_out, body, jd0, T_p, N, M, apsis_flag,
%             mu_sun, search_nu, vinfi_req, vinfo_req, res_flag)
%   rp_GA1: periapsis radius at GA1 [km]
%   rp_GA2: periapsis radius at GA2 [km]
%
% References:
%   [-] n/a
%
% See also:
%   evaluateTotalResonantCost, wrapperResonantCostForFminsearch, findOptimalDSMParameters,
%   grid_scan_3d, scan_and_refine_1d
%
% Adria Sola Foixench
% April 2026

    % Default options
    if isempty(opts)
        opts = struct();
    end

    % Planetary parameters
    [mu_planet, vmr] = GetBodyProps(body);
    vmr_safety = 1.05 * vmr; % Safety margin on minimum flyby radius

    % Orbital period of the resonant body
    [sma_planet, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(body, jd0);
    T_p = 2*pi*sqrt(sma_planet^3 / mu_sun);


    % Search bounds for vinf_out based on available v-infinity magnitudes
    vinfin = norm(vinfi_req);
    vinfon_required = norm(vinfo_req);

    % Anti-trivial floors: 0.5 km/s absolute, 1.0 km/s for 1:1 resonances
    % N==M is the only case where the co-orbital trivial vinf~=0 solution
    % is geometrically admissible and must be excluded explicitly
    if N == M
        v_floor_abs = 1.0;
        v_floor_frac = 0.3;
    else
        v_floor_abs = 0.5;
        v_floor_frac = 0.3;
    end

    if res_flag == 1
        vmag_min = max(v_floor_abs, vinfon_required * v_floor_frac);
        vmag_max = max(5.0, vinfon_required * 1.5);
        v_ref = vinfo_req;
    elseif res_flag == 3
        vmag_min = max(v_floor_abs, vinfin * v_floor_frac);
        vmag_max = max(5.0, vinfin * 1.5);
        v_ref = vinfi_req;
    else
        vmag_min = max(v_floor_abs, min(vinfin, vinfon_required) * v_floor_frac);
        vmag_max = max(5.0, max(vinfin, vinfon_required) * 1.5);
        v_ref = vinfi_req;
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
    vmag_vec = linspace(vmag_min, vmag_max, opts.N_mag);
    theta_vec = linspace(theta_ref - pi, theta_ref + pi, opts.N_theta + 1);
    theta_vec = theta_vec(1:end-1); % Remove duplicate endpoint

    if in_plane
        phi_vec = 0;
    else
        phi_vec = linspace(phi_ref - pi/4, phi_ref + pi/4, opts.N_phi);
    end

    % Define the 1D nu search method
    search_nu = @(f, a, b) scan_and_refine_1d(f, a, b, opts.N_points_nu, opts.N_refine_nu);

    % Genethic Algorithm 
    if strcmpi(opts.algorithm,'ga') || strcmpi(opts.algorithm,'ga+fmin')
        scn.body = body; scn.jd0 = jd0; scn.N = N; scn.M = M; scn.apsis_flag = apsis_flag; scn.res_flag = res_flag; scn.vinfi_req = vinfi_req; scn.vinfo_req = vinfo_req;
        scn.mu_sun = mu_sun; scn.mu_planet = mu_planet; scn.vmr_safety = vmr_safety; scn.T_p = T_p; scn.penaltyCost = 1e6;

        ranges = [vmag_min, vmag_max; rad2deg(theta_ref - pi), rad2deg(theta_ref + pi); rad2deg(phi_ref - pi/4), rad2deg(phi_ref + pi/4); 0, 360; 0, 1];
      
        if in_plane, ranges(3,:) = [0 0]; end

        mutationStep = [0.2, 8, 3, 8];

        fitfun = @(ind) aga_fun_cost(ind, scn, ranges);
        ranfun = @() aga_fun_random(ranges);
        mutfun = @(ind, fit) aga_fun_mutate(ind, fit, ranges, mutationStep);
        repfun = @(a,b,fa,fb) aga_fun_crossover(a, b, fa, fb, ranges);
        unifun = @(pop,fit) aga_fun_unique(pop, fit);

        if isfield(opts,'ga_np')   
            np_ga = opts.ga_np;    
        else
            np_ga = 20; 
        end
        if isfield(opts,'ga_ng')
            ng_ga = opts.ga_ng;    
        else 
            ng_ga = 15; 
        end
        if isfield(opts,'ga_Nctrl')
            Nctrl = opts.ga_Nctrl; 
        else
            Nctrl = [4, floor(np_ga*0.25), floor(np_ga*0.10), floor(np_ga*0.30),]; 
        end

        opts_aga = struct('ninfo',0,'label',0,'dopar',0,'nhist',1);

        % Smart seeding
        n_top = max(2, round(0.15 * np_ga));
        n_jit = max(2, round(0.15 * np_ga));
        max_informed = floor(0.40 * np_ga);

        if n_top + n_jit > max_informed
            n_jit = max(0, max_informed - n_top);
        end

        vmag_c = linspace(vmag_min, vmag_max, 3);
        theta_c = linspace(rad2deg(theta_ref)-180, rad2deg(theta_ref)+180, 6);
        if in_plane
            phi_c = 0;
        else
            phi_c = [rad2deg(phi_ref-pi/4), rad2deg(phi_ref+pi/4)];
        end
        nu_c = [0, 180];

        candidates = {};
        for vv = vmag_c
            for tt = theta_c
                for pp = phi_c
                    for nn = nu_c
                        candidates{end+1} = [vv, tt, pp, nn, 0]; 
                    end
                end
            end
        end
        costs = cellfun(@(c) aga_fun_cost(c, scn, ranges), candidates);

        seeds = cell(1, numel(nu_c));
        for j = 1:numel(nu_c)
            same_nu = cellfun(@(c) c(4) == nu_c(j), candidates);
            sub_idx = find(same_nu);
            [~, b] = min(costs(same_nu));
            seeds{j} = candidates{sub_idx(b)};
        end
        n_top = numel(seeds);

        if n_top + n_jit > max_informed
            n_jit = max(0, max_informed - n_top);
        end

        jitter = [0.1, 2, 2, 2, 0];
        pop0 = cell(1, np_ga);
        k = 1;
        for i = 1:n_top
            if k > np_ga, break; end
            pop0{k} = seeds{i}; k = k + 1;
        end
        for i = 1:n_jit
            if k > np_ga, break; end
            parent = seeds{mod(i-1, n_top) + 1};
            child = parent + jitter .* randn(1, 5);
            child(5) = double(rand >= 0.5);
            pop0{k} = child; k = k + 1;
        end
        while k <= np_ga
            pop0{k} = ranfun(); k = k + 1;
        end
        init_pop = pop0;

        [bestind, ~] = aga(opts_aga, init_pop, 1e-5, ng_ga, Nctrl, unifun, fitfun, mutfun, repfun, ranfun, []);

        % Local Refinment
        if strcmpi(opts.algorithm,'ga+fmin')
            lp_fix = round(bestind(5));
            obj = @(x) aga_fun_cost([x, lp_fix], scn, ranges);
            options_fmin = optimset('TolX',1e-6,'TolFun',1e-6, ...
                'MaxIter',500,'Display','off');
            x_opt = fminsearch(obj, bestind(1:4), options_fmin);
            bestind = [x_opt, lp_fix];
        end

        % Reconstruct breakdown from the GA chromosome
        [~, dec] = aga_fun_decode(bestind, ranges);
        best_vinf_out = dec.vinf_out;
        [~, dbg] = aga_fun_cost(bestind, scn, ranges);
        dV_total = dbg.totalDv;
        best_dV_GA1 = dbg.dV_GA1;
        best_dV_DSM = dbg.dV_DSM;
        best_dV_GA2 = dbg.dV_GA2;
        best_vinf_in = dbg.vinf_in;
        best_va = dbg.va;
        rp_GA1 = dbg.rp1;
        rp_GA2 = dbg.rp2;
        best_vinf_o_n = norm(best_vinf_out);

        vilm_arc = struct('vinf_out',best_vinf_out,'body',body,'jd0',jd0, 'T_p',T_p,'N',N,'M',M,'apsis_flag',apsis_flag,'mu_sun',mu_sun, 'search_nu',search_nu,'vinfi_req',vinfi_req, 'vinfo_req',vinfo_req,'res_flag',res_flag);
        return;
    end

    % 3D grid scan
    f_scan = @(x) evaluateTotalResonantCost(x(1) * [cos(x(3))*cos(x(2)), cos(x(3))*sin(x(2)), sin(x(3))], vinfi_req, vinfo_req, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, search_nu, false);

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
        vilm_arc = struct('vinf_out',[NaN NaN NaN],'body',body,'jd0',jd0, 'T_p',T_p,'N',N,'M',M,'apsis_flag',apsis_flag,'mu_sun',mu_sun, 'search_nu',[],'vinfi_req',vinfi_req,'vinfo_req',vinfo_req,'res_flag',res_flag);
        rp_GA1 = NaN;
        rp_GA2 = NaN;
        return;
    end


    refine = strcmpi(opts.algorithm, 'grid+fmin');

    if refine
        % Local refinement with fminsearch starting from the scan optimum
        [theta0, phi0, vmag0] = cart2sph(vinf_out_scan(1), vinf_out_scan(2), vinf_out_scan(3));

        options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 500, 'Display', 'off');

        if in_plane
            % 2D: optimize only [vmag, theta]
            obj_fmin = @(x) wrapperResonantCostForFminsearch([x(1), x(2), 0], vinfi_req, vinfo_req, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, vmag_min, vmag_max, res_flag, search_nu);
            x_opt2d = fminsearch(obj_fmin, [vmag0, theta0], options);
            x_opt = [x_opt2d(1), x_opt2d(2), 0];
        else
            obj_fmin = @(x) wrapperResonantCostForFminsearch(x, vinfi_req, vinfo_req, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, vmag_min, vmag_max, res_flag, search_nu);
            x0 = [vmag0, theta0, phi0];
            x_opt = fminsearch(obj_fmin, x0, options);
        end

        % Reconstruct refined cartesian vector
        vmag_opt = max(vmag_min, min(vmag_max, x_opt(1)));
        theta_opt = x_opt(2);
        phi_opt = x_opt(3);
        vinf_out_refined = vmag_opt * [cos(phi_opt)*cos(theta_opt), cos(phi_opt)*sin(theta_opt), sin(phi_opt)];

        % Compare refined vs scan and select best
        [dV_refined, dV_GA1_r, dV_DSM_r, dV_GA2_r, vinf_in_r, va_r, ~, ~, ~, rp1_r, rp2_r] = evaluateTotalResonantCost(vinf_out_refined, vinfi_req, vinfo_req, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, search_nu, true);

        if dV_refined < dV_best_scan
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
            [dV_scan_full, dV_GA1_s, dV_DSM_s, dV_GA2_s, vinf_in_s, va_s, ~, ~, ~, rp1_s, rp2_s] = evaluateTotalResonantCost(vinf_out_scan, vinfi_req, vinfo_req, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, search_nu, true);

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
    else
        % Grid only: re-evaluate scan optimum in full mode to recover the
        % full breakdown and trajectory vectors
        [dV_scan_full, dV_GA1_s, dV_DSM_s, dV_GA2_s, vinf_in_s, va_s, ~, ~, ~, rp1_s, rp2_s] = evaluateTotalResonantCost(vinf_out_scan, vinfi_req, vinfo_req, body, jd0, T_p, N, M, apsis_flag, mu_sun, mu_planet, vmr_safety, res_flag, search_nu, true);

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

    best_vinf_o_n = norm(best_vinf_out);

    % Pack inputs needed to rebuild the resonant arcs and to replay the
    % optimizer in validation tests
    vilm_arc.vinf_out = best_vinf_out;
    vilm_arc.body = body;
    vilm_arc.jd0 = jd0;
    vilm_arc.T_p = T_p;
    vilm_arc.N = N;
    vilm_arc.M = M;
    vilm_arc.apsis_flag = apsis_flag;
    vilm_arc.mu_sun = mu_sun;
    vilm_arc.search_nu = search_nu;
    vilm_arc.vinfi_req = vinfi_req;
    vilm_arc.vinfo_req = vinfo_req;
    vilm_arc.res_flag = res_flag;

end
