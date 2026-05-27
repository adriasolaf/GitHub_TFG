function [outerMap, bestRow] = runOuterSearchSpaceMap(config, options)
%runOuterSearchSpaceMap Outer resonant search map
%   Sweeps outgoing v-infinity magnitude and direction. Each sampled point
%   delegates DSM placement to the resonant inner optimizer.
%
% Inputs:
%   config: normalized ROGUI mission/search configuration structure
%   options: optional callbacks structure with cancelFcn and progressFcn
%
% Outputs:
%   outerMap: table with sampled outer candidates and cost diagnostics
%   bestRow: best feasible row from outerMap
%
% Example:
%   [ outerMap, bestRow ] = runOuterSearchSpaceMap ( config );
%
% References:
%   [-]
%
%May 2026

    if nargin < 2
        options = struct();
    end
    options = applyOuterDefaults(options);

    vmagVec = config.outerVmag(:)';
    thetaVec = config.outerTheta(:)';
    phiVec = config.outerPhi(:)';
    nRows = numel(vmagVec) * numel(thetaVec) * numel(phiVec);

    % Preallocate a struct array because table row growth is expensive in
    % MATLAB. Convert to table only after the sweep finishes or is canceled.
    rows = repmat(emptyOuterRow(), nRows, 1);
    rowIdx = 0;

    for vmag = vmagVec
        for theta = thetaVec
            for phi = phiVec
                rowIdx = rowIdx + 1;
                if callCancel(options)
                    rows = rows(1:rowIdx - 1);
                    outerMap = struct2table(rows);
                    bestRow = pickBestOuter(outerMap);
                    return;
                end

                % The outer search variables are spherical coordinates of
                % the planet-relative outgoing excess velocity.
                vinf_out = vmag * [cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi)];
                rows(rowIdx) = evaluateOuterPoint(config, vmag, theta, phi, vinf_out);
                callProgress(options, rowIdx, nRows, rows(rowIdx));
            end
        end
    end

    outerMap = struct2table(rows);
    bestRow = pickBestOuter(outerMap);
end

function row = evaluateOuterPoint(config, vmag, theta, phi, vinf_out)
    row = emptyOuterRow();
    row.vmag = vmag;
    row.theta = theta;
    row.phi = phi;
    row.vinfX = vinf_out(1);
    row.vinfY = vinf_out(2);
    row.vinfZ = vinf_out(3);

    try
        % Auto branch mode uses the repository's reference optimizer. Fixed
        % lp mode uses the local equivalent below so the GUI can enforce a
        % selected Lambert period branch without changing core functions.
        if isempty(config.fixedLp)
            [dV_total, dV_GA1, dV_DSM, dV_GA2, vinf_in, va, r_m, ~, ~, rp1, rp2] = evaluateTotalResonantCost( ...
                vinf_out, config.vinfi, config.vinfo, config.body, config.jd0, config.T_p, config.n_val, config.m_val, ...
                config.apsis_flag, config.mu_sun, config.mu_planet, config.vmr_safety, config.res_flag, config.search_nu, true);

            [dV_best, nu_best, revs_best, lp_best] = findOptimalDSMParameters( ...
                vinf_out, config.body, config.jd0, config.T_p, config.n_val, config.m_val, ...
                config.apsis_flag, config.mu_sun, config.search_nu, false);
        else
            [dV_total, dV_GA1, dV_DSM, dV_GA2, vinf_in, va, r_m, ~, ~, rp1, rp2, ...
                dV_best, nu_best, revs_best, lp_best] = evaluateFixedBranchOuterPoint(config, vinf_out);
        end

        if ~isfinite(dV_total)
            row.failureReason = "No feasible inner optimum.";
            return;
        end

        row.totalDv = dV_total;
        row.dV_GA1 = dV_GA1;
        row.dV_DSM = dV_DSM;
        row.dV_GA2 = dV_GA2;
        row.nuBest = nu_best;
        row.revsBest = revs_best;
        row.lpBest = lp_best;
        row.rp1 = rp1;
        row.rp2 = rp2;
        row.vinfInX = vinf_in(1);
        row.vinfInY = vinf_in(2);
        row.vinfInZ = vinf_in(3);
        row.vaX = va(1);
        row.vaY = va(2);
        row.vaZ = va(3);
        row.rDsmX = r_m(1);
        row.rDsmY = r_m(2);
        row.rDsmZ = r_m(3);
        row.isFeasible = isfinite(dV_best);
        row.failureReason = "";
        row.cacheKey = string(buildSearchSpaceCacheKey(config, row));
    catch err
        row.failureReason = string(err.message);
    end
end

function [dV_total, dV_GA1, dV_DSM, dV_GA2, vinf_in, va, r_m, v_m_minus, v_m_plus, rp1, rp2, ...
    dV_best, nu_best, revs_best, lp_best] = evaluateFixedBranchOuterPoint(config, vinf_out)
    % Local clone of the inner optimizer with branch filters. This keeps
    % fixed-lp behavior in ROGUI and avoids adding GUI concerns to the
    % scientific solver functions.

    dV_total = Inf;
    dV_GA1 = NaN;
    dV_DSM = NaN;
    dV_GA2 = NaN;
    vinf_in = [NaN NaN NaN];
    va = [NaN NaN NaN];
    r_m = [NaN NaN NaN];
    v_m_minus = [NaN NaN NaN];
    v_m_plus = [NaN NaN NaN];
    rp1 = NaN;
    rp2 = NaN;
    dV_best = Inf;
    nu_best = NaN;
    revs_best = NaN;
    lp_best = NaN;

    [model, status] = prepareInnerDsmModel(config, vinf_out);
    if ~status.ok
        return;
    end

    % Enumerate the same revolution splits and Lambert branches as
    % findOptimalDSMParameters, then narrow lp if the GUI fixed it.
    for revs_try = (config.n_val - 1):-1:0
        mr_lambert = config.n_val - revs_try - 1;
        lpValues = branchValues(config.fixedLp, allowedLpValues(mr_lambert));
        for lp_try = lpValues
            % For each discrete branch, minimize only over nu_DSM.
            f_nu = @(nu) computeSingleDSMTransferCost(nu, revs_try, lp_try, config.n_val, ...
                config.apsis_flag, config.mu_sun, false, model.orb_init, model.planets_state);
            nu_try = config.search_nu(f_nu, 0, pi - 1e-3);
            if isnan(nu_try)
                continue;
            end
            dV_try = computeSingleDSMTransferCost(nu_try, revs_try, lp_try, config.n_val, ...
                config.apsis_flag, config.mu_sun, false, model.orb_init, model.planets_state);
            if dV_try < dV_best
                dV_best = dV_try;
                nu_best = nu_try;
                revs_best = revs_try;
                lp_best = lp_try;
            end
        end
    end

    if ~isfinite(dV_best)
        return;
    end

    % Recompute the best branch in full-output mode so the GUI can plot the
    % DSM point, arrival state, and selected trajectory.
    [dV_DSM, r_m, v_m_minus, v_m_plus, vinf_in, va] = computeSingleDSMTransferCost(nu_best, revs_best, lp_best, ...
        config.n_val, config.apsis_flag, config.mu_sun, true, model.orb_init, model.planets_state);
    [dV_GA1, rp1, ok1] = gravityAssistCost(config.vinfi, vinf_out, config.mu_planet, config.vmr_safety, config.res_flag == 1);
    [dV_GA2, rp2, ok2] = gravityAssistCost(vinf_in, config.vinfo, config.mu_planet, config.vmr_safety, config.res_flag == 3);
    if ok1 && ok2
        dV_total = dV_GA1 + dV_DSM + dV_GA2;
    end
end

function values = branchValues(fixedValue, autoValues)
    if isempty(fixedValue)
        values = autoValues;
    elseif any(autoValues == fixedValue)
        values = fixedValue;
    else
        values = [];
    end
end

function values = allowedLpValues(mr_lambert)
    if mr_lambert > 0
        values = [0, 1];
    else
        values = 0;
    end
end

function [dV, rp, ok] = gravityAssistCost(vinf_in, vinf_out, mu_planet, vmr_safety, skipCost)
    % Same patched-conics powered flyby cost used by the core code, with
    % cosine clamping to avoid numerical acos drift outside [-1, 1].
    dV = 0;
    rp = NaN;
    ok = true;
    if skipCost
        return;
    end
    vinfin = norm(vinf_in);
    vinfon = norm(vinf_out);
    if vinfin < 1e-9 || vinfon < 1e-9
        ok = false;
        dV = NaN;
        return;
    end
    cosDelta = dot(vinf_in, vinf_out) / (vinfin * vinfon);
    cosDelta = max(-1.0, min(1.0, cosDelta));
    delta = acos(cosDelta);
    [dV, rp] = GA_PGA2_Rp(vinfin, vinfon, delta, mu_planet);
    if rp < vmr_safety
        rp = vmr_safety;
        dV = GA_PGA2_Vinfo(vinfin, vinfon, delta, vmr_safety, mu_planet);
    end
end

function row = emptyOuterRow()
    row = struct('vmag', NaN, 'theta', NaN, 'phi', NaN, ...
        'vinfX', NaN, 'vinfY', NaN, 'vinfZ', NaN, ...
        'totalDv', Inf, 'dV_GA1', NaN, 'dV_DSM', NaN, 'dV_GA2', NaN, ...
        'nuBest', NaN, 'revsBest', NaN, 'lpBest', NaN, ...
        'rp1', NaN, 'rp2', NaN, ...
        'vinfInX', NaN, 'vinfInY', NaN, 'vinfInZ', NaN, ...
        'vaX', NaN, 'vaY', NaN, 'vaZ', NaN, ...
        'rDsmX', NaN, 'rDsmY', NaN, 'rDsmZ', NaN, ...
        'isFeasible', false, 'failureReason', "Not evaluated.", 'cacheKey', "");
end

function bestRow = pickBestOuter(outerMap)
    bestRow = table();
    if isempty(outerMap)
        return;
    end
    feasible = outerMap.isFeasible & isfinite(outerMap.totalDv);
    if ~any(feasible)
        return;
    end
    feasibleRows = outerMap(feasible, :);
    [~, idx] = min(feasibleRows.totalDv);
    bestRow = feasibleRows(idx, :);
end

function options = applyOuterDefaults(options)
    if ~isfield(options, 'cancelFcn'), options.cancelFcn = []; end
    if ~isfield(options, 'progressFcn'), options.progressFcn = []; end
end

function tf = callCancel(options)
    tf = false;
    if ~isempty(options.cancelFcn)
        tf = logical(options.cancelFcn());
    end
end

function callProgress(options, idx, total, row)
    if ~isempty(options.progressFcn)
        options.progressFcn(idx, total, row);
    end
    if mod(idx, 10) == 0
        drawnow limitrate;
    end
end
