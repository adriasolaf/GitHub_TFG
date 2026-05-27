function [ cost, dbg ] = aga_fun_cost ( ind, scn, ranges )
%aga_fun_cost Resonant branch AGA cost
%   Evaluates one fixed [vinfo,theta,phi,nu,lp] chromosome. The best
%   feasible revolution split is selected internally for the fixed nu/lp.
%
% Inputs:
%   ind: AGA individual [vinfo_km_s,theta_deg,phi_deg,nu_deg,lp]
%   scn: resonant branch scenario structure with fields:
%       body, jd0, N, M, apsis_flag, res_flag, vinfi_req, vinfo_req
%       mu_sun, mu_planet, vmr_safety, T_p
%   ranges: chromosome bounds, one [min,max] row per gene
%
% Outputs:
%   cost: total resonant branch cost [km s^-1]
%   dbg: debug data structure
%
% Example:
%   [ cost, dbg ] = aga_fun_cost ( ind, scn, ranges );
%
% References:
%   [-]
%
%May 2026

scn = normalizeScenario(scn);
penaltyCost = getOption(scn, 'penaltyCost', 1e12);
if nargin < 3 || isempty(ranges)
    ranges = aga_aux_default_ranges();
end

% Decode once. The decoder clamps genes to ranges, converts angles to
% radians, and builds the Cartesian outgoing v-infinity vector.
[~, decoded] = aga_fun_decode(ind, ranges);
dbg = aga_aux_empty_dbg(decoded);
cost = penaltyCost;
if ~decoded.isValid || decoded.vinfo < 1e-9
    dbg.failureReason = "Invalid GA individual.";
    return;
end

try
    % Build all Kepler/Lambert inputs that depend on this chromosome's
    % outgoing excess velocity. Invalid heliocentric departure orbits are
    % caught and converted into the penalty value below.
    model = buildDsmModel(scn, decoded.vinf_out);

    % Keep nu/lp fixed from the chromosome. The current core evaluator also
    % needs revs_before, so scan the admissible splits and keep the best.
    [dV_DSM, revs_best, r_m, v_m_minus, v_m_plus, vinf_in, va] = evaluateFixedNuLp(scn, model, decoded.nu, decoded.lp);
    if ~isfinite(dV_DSM)
        dbg.failureReason = "No feasible DSM/Lambert solution for fixed nu/lp.";
        return;
    end

    % Add the powered flyby terms adjacent to the resonant leg. res_flag
    % skips the missing edge flyby for launch-side or arrival-side cases.
    [dV_GA1, rp1, ok1] = gravityAssistCost(getScenarioVector(scn, 'vinfi_req', 'vinfi'), decoded.vinf_out, ...
        scn.mu_planet, scn.vmr_safety, getOption(scn, 'res_flag', 2) == 1);
    [dV_GA2, rp2, ok2] = gravityAssistCost(vinf_in, getScenarioVector(scn, 'vinfo_req', 'vinfo'), ...
        scn.mu_planet, scn.vmr_safety, getOption(scn, 'res_flag', 2) == 3);
    if ~(ok1 && ok2)
        dbg.failureReason = "Invalid powered flyby geometry.";
        return;
    end

    % Return the scalar objective used by AGA and keep useful debug data for
    % post-run inspection.
    cost = dV_GA1 + dV_DSM + dV_GA2;
    dbg.isFeasible = true;
    dbg.failureReason = "";
    dbg.totalDv = cost;
    dbg.dV_GA1 = dV_GA1;
    dbg.dV_DSM = dV_DSM;
    dbg.dV_GA2 = dV_GA2;
    dbg.revs_before = revs_best;
    dbg.rp1 = rp1;
    dbg.rp2 = rp2;
    dbg.r_m = r_m;
    dbg.v_m_minus = v_m_minus;
    dbg.v_m_plus = v_m_plus;
    dbg.vinf_in = vinf_in;
    dbg.va = va;
catch err
    dbg.failureReason = string(err.message);
end

end

function [ scn ] = normalizeScenario ( scn )
%normalizeScenario Complete scenario aliases and constants.

% Accept both core names and ROGUI-derived names. This lets the same cost
% function work with a hand-built scenario and with a GUI-derived struct.
if ~isfield(scn, 'N') && isfield(scn, 'n_val')
    scn.N = scn.n_val;
end
if ~isfield(scn, 'M') && isfield(scn, 'm_val')
    scn.M = scn.m_val;
end
if ~isfield(scn, 'vinfi_req') && isfield(scn, 'vinfi')
    scn.vinfi_req = scn.vinfi;
end
if ~isfield(scn, 'vinfo_req') && isfield(scn, 'vinfo')
    scn.vinfo_req = scn.vinfo;
end

% Fill physical constants when the caller provides a lightweight scenario.
% The explicit example script already fills these, so this is only fallback.
if ~isfield(scn, 'mu_sun') || isempty(scn.mu_sun)
    scn.mu_sun = GetBodyProps('Sun');
end
if ~isfield(scn, 'mu_planet') || ~isfield(scn, 'vmr_safety')
    [muPlanet, radius] = GetBodyProps(scn.body);
    if ~isfield(scn, 'mu_planet') || isempty(scn.mu_planet)
        scn.mu_planet = muPlanet;
    end
    if ~isfield(scn, 'vmr_safety') || isempty(scn.vmr_safety)
        scn.vmr_safety = 1.05 * radius;
    end
end
if ~isfield(scn, 'T_p') || isempty(scn.T_p)
    [smaPlanet, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(scn.body, scn.jd0);
    scn.T_p = 2*pi*sqrt(smaPlanet^3 / scn.mu_sun);
end

end

function [ model ] = buildDsmModel ( scn, vinf_out )
%buildDsmModel Build DSM state model for one chromosome.

% The resonant arrival epoch is M revolutions of the resonant body.
sec2days = 1 / 86400;
model.tof_total_s = scn.M * scn.T_p;
model.jd_f = scn.jd0 + model.tof_total_s * sec2days;

[model.r_p0, model.v_p0] = GetBodyICF(scn.body, scn.jd0, scn.mu_sun, 0);
[model.r_pf, model.v_pf] = GetBodyICF(scn.body, model.jd_f, scn.mu_sun, 0);

model.planets_state.r_pf = model.r_pf;
model.planets_state.v_pf = model.v_pf;
model.planets_state.jd_f = model.jd_f;
model.planets_state.tof_total_s = model.tof_total_s;
model.v_sc0 = model.v_p0 + vinf_out;

% The pre-DSM orbit must be elliptic for this resonant model.
[sma, ecc, inc, nu0, argp, raan] = ICF2KEP_O(model.r_p0, model.v_sc0, scn.mu_sun);
if ecc >= 1.0 || sma <= 0
    error('AGA:InvalidOrbit', 'Candidate produces a hyperbolic or invalid heliocentric orbit.');
end

% Store the compact orbit state expected by computeSingleDSMTransferCost.
% M0 anchors the time from departure to the chromosome-selected DSM anomaly.
n_motion = sqrt(scn.mu_sun / sma^3);
M0 = trueAnomalyToMeanAnomaly(nu0, ecc);

model.orb_init.sma = sma;
model.orb_init.ecc = ecc;
model.orb_init.inc = inc;
model.orb_init.argp = argp;
model.orb_init.raan = raan;
model.orb_init.n_motion = n_motion;
model.orb_init.M0 = M0;
model.orb_init.v_sc0 = model.v_sc0;

end

function [dV_best, revs_best, r_m, v_m_minus, v_m_plus, vinf_in, va] = evaluateFixedNuLp ( scn, model, nu, lp )
%evaluateFixedNuLp Evaluate one fixed nu/lp pair.

% Initialize as infeasible. Any valid candidate below will replace these.
dV_best = Inf;
revs_best = NaN;
r_m = [NaN NaN NaN];
v_m_minus = [NaN NaN NaN];
v_m_plus = [NaN NaN NaN];
vinf_in = [NaN NaN NaN];
va = [NaN NaN NaN];

for revs_try = (scn.N - 1):-1:0
    mr_lambert = scn.N - revs_try - 1;

    % Long-period Lambert is only meaningful for multi-revolution arcs.
    if mr_lambert == 0 && lp == 1
        continue;
    end

    % Evaluate the selected DSM anomaly and Lambert branch for this split.
    [dV_try, r_try, vm_minus_try, vm_plus_try, vinf_in_try, va_try] = computeSingleDSMTransferCost( ...
        nu, revs_try, lp, scn.N, scn.apsis_flag, scn.mu_sun, true, model.orb_init, model.planets_state);
    if dV_try < dV_best
        dV_best = dV_try;
        revs_best = revs_try;
        r_m = r_try;
        v_m_minus = vm_minus_try;
        v_m_plus = vm_plus_try;
        vinf_in = vinf_in_try;
        va = va_try;
    end
end

end

function [ dV, rp, ok ] = gravityAssistCost ( vinf_in, vinf_out, mu_planet, vmr_safety, skipCost )
%gravityAssistCost Compute one powered flyby cost.

% Edge resonant branches do not have both adjacent flybys.
dV = 0;
rp = NaN;
ok = true;
if skipCost
    return;
end

vinfin = norm(vinf_in);
vinfon = norm(vinf_out);
if vinfin < 1e-9 || vinfon < 1e-9
    dV = NaN;
    ok = false;
    return;
end

cosDelta = dot(vinf_in, vinf_out) / (vinfin * vinfon);
cosDelta = max(-1.0, min(1.0, cosDelta));
delta = acos(cosDelta);

% First compute the natural powered flyby radius, then enforce safety.
[dV, rp] = GA_PGA2_Rp(vinfin, vinfon, delta, mu_planet);
if rp < vmr_safety
    rp = vmr_safety;
    dV = GA_PGA2_Vinfo(vinfin, vinfon, delta, vmr_safety, mu_planet);
end

end

function [ value ] = getOption ( s, name, defaultValue )
%getOption Return a struct field or a default value.

if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end

end

function [ value ] = getScenarioVector ( scn, primaryName, fallbackName )
%getScenarioVector Return a scenario vector field as a row.

if isfield(scn, primaryName)
    value = scn.(primaryName);
else
    value = scn.(fallbackName);
end
value = value(:)';

end
