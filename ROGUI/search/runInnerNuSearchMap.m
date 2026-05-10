function [innerMap, bestRow] = runInnerNuSearchMap(config, outerCandidate, options)
%runInnerNuSearchMap Evaluate DSM nu values and Lambert branches.
%   This exposes the inner search space for a fixed outgoing vinf vector.
%
%   The returned table is intentionally branch-resolved: the same nu_DSM
%   can appear multiple times with different revs_before/lw/lp values.
%   This is why the GUI can show a double-valued or multi-valued curve.

    if nargin < 3
        options = struct();
    end
    options = applyInnerDefaults(options);

    vinf_out = extractVinfOut(outerCandidate);
    [model, status] = prepareInnerDsmModel(config, vinf_out);
    if ~status.ok
        innerMap = struct2table(emptyInnerRow());
        innerMap.failureReason(1) = string(status.message);
        bestRow = table();
        return;
    end

    % innerNu is a GUI-controlled display grid, not the optimizer's refine
    % grid. It shows what the inner model looks like across sampled nu.
    nuVec = config.innerNu(:)';
    maxRows = numel(nuVec) * config.n_val * 4;
    rows = repmat(emptyInnerRow(), maxRows, 1);
    rowIdx = 0;

    % Match the branch families explored by findOptimalDSMParameters while
    % honoring optional fixed lw/lp values from the GUI.
    for revs_before = (config.n_val - 1):-1:0
        mr_lambert = config.n_val - revs_before - 1;
        lwValues = branchValues(config.fixedLw, [0, 1]);
        lpValues = branchValues(config.fixedLp, allowedLpValues(mr_lambert));

        if isempty(lpValues)
            continue;
        end

        for lw = lwValues
            for lp = lpValues
                for nu = nuVec
                    rowIdx = rowIdx + 1;
                    if callCancel(options)
                        rows = rows(1:rowIdx - 1);
                        innerMap = struct2table(rows);
                        bestRow = pickBestInner(innerMap);
                        return;
                    end
                    rows(rowIdx) = evaluateInnerPoint(config, model, vinf_out, nu, revs_before, lw, lp);
                    callProgress(options, rowIdx, maxRows, rows(rowIdx));
                end
            end
        end
    end

    rows = rows(1:rowIdx);
    innerMap = struct2table(rows);
    bestRow = pickBestInner(innerMap);
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

function row = evaluateInnerPoint(config, model, vinf_out, nu, revs_before, lw, lp)
    % Evaluate one visible sample in the inner plot. This is not an
    % optimizer step; it is a diagnostic evaluation of a specific branch
    % and nu value.
    row = emptyInnerRow();
    row.nu = nu;
    row.revs_before = revs_before;
    row.lw = lw;
    row.lp = lp;

    try
        [dV_DSM, r_m, v_m_minus, v_m_plus, vinf_in, va] = computeSingleDSMTransferCost( ...
            nu, revs_before, lw, lp, config.n_val, config.apsis_flag, config.mu_sun, true, model.orb_init, model.planets_state);

        if ~isfinite(dV_DSM)
            row.failureReason = "Invalid DSM timing or Lambert branch.";
            return;
        end

        [dV_GA1, rp1, ok1] = gravityAssistCost(config.vinfi, vinf_out, config.mu_planet, config.vmr_safety, config.res_flag == 1);
        [dV_GA2, rp2, ok2] = gravityAssistCost(vinf_in, config.vinfo, config.mu_planet, config.vmr_safety, config.res_flag == 3);

        if ~(ok1 && ok2)
            row.failureReason = "Invalid gravity-assist geometry.";
            return;
        end

        row.dV_DSM = dV_DSM;
        row.dV_GA1 = dV_GA1;
        row.dV_GA2 = dV_GA2;
        row.totalDv = dV_GA1 + dV_DSM + dV_GA2;
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
        row.vMinusX = v_m_minus(1);
        row.vMinusY = v_m_minus(2);
        row.vMinusZ = v_m_minus(3);
        row.vPlusX = v_m_plus(1);
        row.vPlusY = v_m_plus(2);
        row.vPlusZ = v_m_plus(3);
        row.isFeasible = true;
        row.failureReason = "";
    catch err
        row.failureReason = string(err.message);
    end
end

function [dV, rp, ok] = gravityAssistCost(vinf_in, vinf_out, mu_planet, vmr_safety, skipCost)
    % Keep the inner-map cost consistent with evaluateTotalResonantCost:
    % DSM plus the powered flyby terms that apply at this resonance position.
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

function vinf_out = extractVinfOut(outerCandidate)
    if istable(outerCandidate)
        outerCandidate = outerCandidate(1, :);
        vinf_out = [outerCandidate.vinfX, outerCandidate.vinfY, outerCandidate.vinfZ];
    elseif isstruct(outerCandidate) && isfield(outerCandidate, 'vinf_out')
        vinf_out = outerCandidate.vinf_out;
    elseif isstruct(outerCandidate) && isfield(outerCandidate, 'vinfX')
        vinf_out = [outerCandidate.vinfX, outerCandidate.vinfY, outerCandidate.vinfZ];
    else
        vinf_out = outerCandidate(:)';
    end
end

function row = emptyInnerRow()
    row = struct('nu', NaN, 'revs_before', NaN, 'lw', NaN, 'lp', NaN, ...
        'totalDv', Inf, 'dV_GA1', NaN, 'dV_DSM', NaN, 'dV_GA2', NaN, ...
        'rp1', NaN, 'rp2', NaN, ...
        'vinfInX', NaN, 'vinfInY', NaN, 'vinfInZ', NaN, ...
        'vaX', NaN, 'vaY', NaN, 'vaZ', NaN, ...
        'rDsmX', NaN, 'rDsmY', NaN, 'rDsmZ', NaN, ...
        'vMinusX', NaN, 'vMinusY', NaN, 'vMinusZ', NaN, ...
        'vPlusX', NaN, 'vPlusY', NaN, 'vPlusZ', NaN, ...
        'isFeasible', false, 'failureReason', "Not evaluated.");
end

function bestRow = pickBestInner(innerMap)
    bestRow = table();
    if isempty(innerMap)
        return;
    end
    feasible = innerMap.isFeasible & isfinite(innerMap.totalDv);
    if ~any(feasible)
        return;
    end
    feasibleRows = innerMap(feasible, :);
    [~, idx] = min(feasibleRows.totalDv);
    bestRow = feasibleRows(idx, :);
end

function options = applyInnerDefaults(options)
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
    if mod(idx, 20) == 0
        drawnow limitrate;
    end
end
