function traj = evaluateTrajectoryCandidate(config, outerCandidate, innerCandidate)
%evaluateTrajectoryCandidate Build trajectory data for selected GUI cases.
%   This reconstructs the configured mission while forcing the selected
%   resonant outgoing vinf and, when supplied, selected DSM branch.
%
%   The returned traj struct mirrors the main MGA2_PGA2 output fields and
%   adds GUI-specific plot data. It is deliberately script-callable so a
%   selected GUI case can be verified outside the interface.

    if nargin < 3
        innerCandidate = [];
    end

    traj = emptyTrajectory();
    vinf_out = extractVinfOut(outerCandidate);

    try
        [inner, innerStatus] = resolveInnerSelection(config, vinf_out, innerCandidate);
        if ~innerStatus.ok
            traj.status.message = innerStatus.message;
            return;
        end

        traj = buildPatchedTrajectory(config, vinf_out, inner);
        traj.status.ok = true;
        traj.status.message = 'OK';
    catch err
        traj.status.ok = false;
        traj.status.message = err.message;
    end
end

function [inner, status] = resolveInnerSelection(config, vinf_out, innerCandidate)
    % If the user has not clicked the inner plot, use the optimizer's best
    % branch. If they have clicked the inner plot, force that exact
    % nu/revs/lw/lp combination for trajectory reconstruction.
    status = struct('ok', false, 'message', '');
    inner = struct();

    [model, modelStatus] = prepareInnerDsmModel(config, vinf_out);
    if ~modelStatus.ok
        status.message = modelStatus.message;
        return;
    end

    if isempty(innerCandidate)
        [dV_DSM, nu, revs_before, lw, lp, r_m, v_m_minus, v_m_plus, vinf_in, va] = findOptimalDSMParameters( ...
            vinf_out, config.body, config.jd0, config.T_p, config.n_val, config.m_val, ...
            config.apsis_flag, config.mu_sun, config.search_nu, true);
    else
        innerCandidate = normalizeCandidate(innerCandidate);
        nu = innerCandidate.nu;
        revs_before = innerCandidate.revs_before;
        lw = innerCandidate.lw;
        lp = innerCandidate.lp;
        [dV_DSM, r_m, v_m_minus, v_m_plus, vinf_in, va] = computeSingleDSMTransferCost( ...
            nu, revs_before, lw, lp, config.n_val, config.apsis_flag, config.mu_sun, true, model.orb_init, model.planets_state);
    end

    if ~isfinite(dV_DSM)
        status.message = 'Selected inner case is infeasible.';
        return;
    end

    inner.dV_DSM = dV_DSM;
    inner.nu = nu;
    inner.revs_before = revs_before;
    inner.lw = lw;
    inner.lp = lp;
    inner.r_m = r_m;
    inner.v_m_minus = v_m_minus;
    inner.v_m_plus = v_m_plus;
    inner.vinf_in = vinf_in;
    inner.va = va;
    inner.model = model;
    status.ok = true;
    status.message = 'OK';
end

function traj = buildPatchedTrajectory(config, vinf_out, inner)
    % Build a patched trajectory without rerunning the outer optimizer.
    % Non-resonant Lambert arcs are already stored in config; the selected
    % resonant leg is injected here using the chosen vinf_out and DSM result.
    traj = emptyTrajectory();
    nPlanets = numel(config.planets);
    nTransfers = nPlanets - 1;
    idxRes = config.resonantIndex;

    jd2k = config.jd2k;
    r = config.r;
    v = config.v;
    vd = config.vd;
    va = config.va;
    dvga = zeros(nTransfers - 1, 1);
    rpga = NaN(nTransfers - 1, 1);
    dvdsm = zeros(nTransfers, 1);

    % Replace only the selected resonant transfer state. Other transfers
    % remain the adjacent Lambert context computed during config creation.
    vd(idxRes, :) = v(idxRes, :) + vinf_out;
    va(idxRes, :) = inner.va;
    dvdsm(idxRes) = inner.dV_DSM;

    [dV_GA1, rp1, ok1] = gravityAssistCost(config.vinfi, vinf_out, config.mu_planet, config.vmr_safety, config.res_flag == 1);
    [dV_GA2, rp2, ok2] = gravityAssistCost(inner.vinf_in, config.vinfo, config.mu_planet, config.vmr_safety, config.res_flag == 3);
    if ~(ok1 && ok2)
        error('ROGUI:InvalidCandidate', 'Selected candidate has invalid gravity-assist geometry.');
    end

    if config.res_flag == 1
        dvga(idxRes) = dV_GA2;
        rpga(idxRes) = rp2;
    elseif config.res_flag == 2
        dvga(idxRes - 1) = dV_GA1;
        rpga(idxRes - 1) = rp1;
        dvga(idxRes) = dV_GA2;
        rpga(idxRes) = rp2;
    else
        dvga(idxRes - 1) = dV_GA1;
        rpga(idxRes - 1) = rp1;
    end

    % Recompute standard powered flyby costs only at non-VILM encounters.
    % VILM-adjacent flybys are accounted for by the selected resonant cost.
    for idx = 1:(nTransfers - 1)
        if config.is_vilm(idx) || config.is_vilm(idx + 1)
            continue;
        end
        [dvga(idx), rpga(idx)] = standardGaCost(va(idx, :) - v(idx + 1, :), vd(idx + 1, :) - v(idx + 1, :), config.planets{idx + 1});
    end

    vilm_arcs = repmat(struct('vinf_out', [], 'body', '', 'jd0', [], 'T_p', [], ...
        'N', [], 'M', [], 'apsis_flag', [], 'mu_sun', [], 'search_nu', []), nTransfers, 1);
    vilm_arcs(idxRes) = struct('vinf_out', vinf_out, 'body', config.body, 'jd0', config.jd0, ...
        'T_p', config.T_p, 'N', config.n_val, 'M', config.m_val, 'apsis_flag', config.apsis_flag, ...
        'mu_sun', config.mu_sun, 'search_nu', config.search_nu);

    traj.jd2k = jd2k;
    traj.r = r;
    traj.v = v;
    traj.vd = vd;
    traj.va = va;
    traj.rpga = rpga;
    traj.dvga = dvga;
    traj.dvdsm = dvdsm;
    traj.vilm_arcs = vilm_arcs;
    traj.planets = config.planets;
    traj.is_vilm = config.is_vilm;
    traj.selectedOuter.vinf_out = vinf_out;
    traj.selectedInner = rmfield(inner, 'model');
    traj.cost.totalDv = dV_GA1 + inner.dV_DSM + dV_GA2;
    traj.cost.dV_GA1 = dV_GA1;
    traj.cost.dV_DSM = inner.dV_DSM;
    traj.cost.dV_GA2 = dV_GA2;
    traj.cost.rp1 = rp1;
    traj.cost.rp2 = rp2;
    traj.plot = buildPlotData(config, traj, inner);
end

function plotData = buildPlotData(config, traj, inner)
    % Convert numerical trajectory outputs into dense point arrays for the
    % plot helper. This keeps plotTrajectory3D free of mission mechanics.
    nPlanets = numel(config.planets);
    nTransfers = nPlanets - 1;
    nOrbit = 150;

    plotData.planetOrbits = cell(nPlanets, 1);
    for idx = 1:nPlanets
        plotData.planetOrbits{idx} = GetBodyOrbit(config.planets{idx}, traj.jd2k(idx), config.mu_sun, nOrbit, 1);
    end

    plotData.transferArcs = cell(nTransfers, 1);
    for idx = 1:nTransfers
        if config.is_vilm(idx)
            continue;
        end
        plotData.transferArcs{idx} = ICF2Arc(traj.r(idx, :), traj.vd(idx, :), traj.r(idx + 1, :), traj.va(idx, :), config.mu_sun, nOrbit);
    end

    idxRes = config.resonantIndex;
    mr_lambert = config.n_val - inner.revs_before - 1;
    [arc1, arc2, rdsm] = generateResonantTrajectoryPoints(inner.model.r_p0, inner.model.v_sc0, inner.r_m, ...
        inner.v_m_minus, inner.v_m_plus, inner.model.r_pf, inner.va, inner.revs_before, mr_lambert, config.mu_sun);
    plotData.resonant(idxRes).arc1 = arc1;
    plotData.resonant(idxRes).arc2 = arc2;
    plotData.resonant(idxRes).rdsm = rdsm;
end

function [dV, rp, ok] = standardGaCost(vinfi, vinfo, body)
    [mu_planet, vmr] = GetBodyProps(body);
    vmr_safety = 1.05 * vmr;
    [dV, rp, ok] = gravityAssistCost(vinfi, vinfo, mu_planet, vmr_safety, false);
end

function [dV, rp, ok] = gravityAssistCost(vinfi, vinfo, mu_planet, vmr_safety, skipCost)
    % Local copy of the powered flyby cost used for selected trajectory
    % metrics. The skipCost flag handles launch/arrival-edge resonances.
    dV = 0;
    rp = NaN;
    ok = true;
    if skipCost
        return;
    end
    vinfin = norm(vinfi);
    vinfon = norm(vinfo);
    if vinfin < 1e-9 || vinfon < 1e-9
        ok = false;
        dV = NaN;
        return;
    end
    cos_delta = dot(vinfi, vinfo) / (vinfin * vinfon);
    cos_delta = max(-1.0, min(1.0, cos_delta));
    delta = acos(cos_delta);
    [dV, rp] = GA_PGA2_Rp(vinfin, vinfon, delta, mu_planet);
    if rp < vmr_safety
        rp = vmr_safety;
        dV = GA_PGA2_Vinfo(vinfin, vinfon, delta, vmr_safety, mu_planet);
    end
end

function vinf_out = extractVinfOut(candidate)
    if istable(candidate)
        candidate = candidate(1, :);
        vinf_out = [candidate.vinfX, candidate.vinfY, candidate.vinfZ];
    elseif isstruct(candidate) && isfield(candidate, 'vinf_out')
        vinf_out = candidate.vinf_out;
    elseif isstruct(candidate) && isfield(candidate, 'vinfX')
        vinf_out = [candidate.vinfX, candidate.vinfY, candidate.vinfZ];
    else
        vinf_out = candidate(:)';
    end
end

function candidate = normalizeCandidate(candidate)
    if istable(candidate)
        candidate = table2struct(candidate(1, :));
    end
end

function traj = emptyTrajectory()
    traj = struct();
    traj.status = struct('ok', false, 'message', 'Not evaluated.');
    traj.jd2k = [];
    traj.r = [];
    traj.v = [];
    traj.vd = [];
    traj.va = [];
    traj.rpga = [];
    traj.dvga = [];
    traj.dvdsm = [];
    traj.vilm_arcs = [];
    traj.planets = {};
    traj.is_vilm = [];
    traj.selectedOuter = struct();
    traj.selectedInner = struct();
    traj.cost = struct();
    traj.plot = struct();
end
