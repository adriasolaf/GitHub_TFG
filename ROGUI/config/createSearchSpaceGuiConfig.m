function config = createSearchSpaceGuiConfig(inputConfig)
%createSearchSpaceGuiConfig Create ROGUI configuration
%   Validates user-facing mission/search inputs and derives the resonant
%   branch data used by ROGUI search and plotting helpers.
%
% Inputs:
%   inputConfig: partial configuration structure
%
% Outputs:
%   config: normalized ROGUI mission/search configuration structure
%
% Example:
%   [ config ] = createSearchSpaceGuiConfig ( inputConfig );
%
% References:
%   [-]
%
%May 2026

    if nargin < 1 || isempty(inputConfig)
        inputConfig = defaultConfig();
    end

    % Start from a complete deterministic Galileo-style case, then replace
    % only the fields provided by the GUI or validation script.
    config = mergeDefaults(inputConfig, defaultConfig());

    if ischar(config.planets) || isstring(config.planets)
        config.planets = parsePlanetSequence(config.planets);
    end
    config.planets = config.planets(:)';

    nPlanets = numel(config.planets);
    if nPlanets < 2
        error('ROGUI:InvalidConfig', 'At least two planets are required.');
    end
    if numel(config.tofs) ~= nPlanets - 1
        error('ROGUI:InvalidConfig', 'TOF vector must contain one value per transfer.');
    end

    config.tofs = config.tofs(:)';
    config.N = normalizeResonanceVector(config.N);
    config.M = normalizeResonanceVector(config.M);

    % The search maps repeatedly evaluate the selected resonant leg. Cache
    % encounter epochs, planet states, and adjacent non-resonant Lambert
    % arcs once in config so downstream functions stay simple.
    config.mu_sun = GetBodyProps('Sun');
    config.days2secs = 86400;
    config.jd2k = zeros(nPlanets, 1);
    config.r = zeros(nPlanets, 3);
    config.v = zeros(nPlanets, 3);
    config.vd = NaN(nPlanets - 1, 3);
    config.va = NaN(nPlanets - 1, 3);

    config.jd2k(1) = config.jd2k0;
    for idx = 2:nPlanets
        config.jd2k(idx) = config.jd2k(idx - 1) + config.tofs(idx - 1);
    end

    for idx = 1:nPlanets
        [config.r(idx, :), config.v(idx, :)] = GetBodyICF(config.planets{idx}, config.jd2k(idx), config.mu_sun, 1);
    end

    [config.N_per_leg, flagN] = scalar2LengthN(config.N, nPlanets - 1);
    [config.M_per_leg, flagM] = scalar2LengthN(config.M, nPlanets - 1);
    if flagN ~= 0 || flagM ~= 0
        error('ROGUI:InvalidConfig', 'N and M must be scalars or vectors with one value per transfer leg.');
    end

    % ROGUI follows the core repository convention: a resonant leg is a
    % same-body transfer whose TOF matches M planetary periods.
    config.is_vilm = false(nPlanets - 1, 1);
    for idx = 1:(nPlanets - 1)
        if ~strcmp(config.planets{idx}, config.planets{idx + 1})
            continue;
        end

        m_val_try = config.M_per_leg(idx);
        if isnan(m_val_try)
            continue;
        end

        [sma_p, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(config.planets{idx}, config.jd2k(idx));
        T_p_s = 2*pi * sqrt(sma_p^3 / config.mu_sun);
        tof_res_days = m_val_try * T_p_s / config.days2secs;
        tol_days = 2;

        config.is_vilm(idx) = abs(config.tofs(idx) - tof_res_days) < tol_days;
    end

    resonantIndices = find(config.is_vilm);
    if isempty(resonantIndices)
        error('ROGUI:InvalidConfig', 'No time-consistent resonant leg found. Use consecutive equal planet names and TOF ~= M*T_p.');
    end
    if isempty(config.resonantIndex)
        config.resonantIndex = resonantIndices(1);
    end
    if ~ismember(config.resonantIndex, resonantIndices)
        error('ROGUI:InvalidConfig', 'Selected resonantIndex is not a time-consistent same-planet resonant transfer.');
    end

    % Scalar N/M values apply to all transfer legs; vector N/M values are
    % indexed by transfer leg, matching MGA2_PGA2.
    config.resonantIndices = resonantIndices;
    config.resonanceOrdinal = find(resonantIndices == config.resonantIndex, 1);
    config.n_val = config.N_per_leg(config.resonantIndex);
    config.m_val = config.M_per_leg(config.resonantIndex);
    config.apsis_flag = inferApsisFlag(config.n_val, config.m_val);
    config.res_flag = inferResFlag(config.resonantIndex, nPlanets - 1);
    config.body = config.planets{config.resonantIndex};
    config.jd0 = config.jd2k(config.resonantIndex);

    [config.mu_planet, config.body_radius] = GetBodyProps(config.body);
    config.vmr_safety = config.flybySafetyFactor * config.body_radius;
    [sma_planet, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(config.body, config.jd0);
    config.T_p = 2*pi*sqrt(sma_planet^3 / config.mu_sun);

    % Store the search settings as data but rebuild the function handle in
    % memory. Function handles are convenient during a session but should
    % not be treated as durable serialized state.
    config.search_nu = @(f, a, b) scan_and_refine_1d(f, a, b, config.innerScanPoints, config.innerRefineIterations);

    config = computeAdjacentLambertContext(config);
end

function config = defaultConfig()
    % Default scenario mirrors TEST_MGA_Plot.m and gives developers a
    % known-good case immediately after opening the GUI.
    config.planets = {'Earth', 'Venus', 'Earth', 'Earth', 'Jupiter'};
    config.jd2k0 = -3727;
    config.tofs = defaultGalileoTofs(config.planets, config.jd2k0, 1, 2);
    config.N = 1;
    config.M = 2;
    config.resonantIndex = [];
    config.flybySafetyFactor = 1.05;
    config.innerScanPoints = 80;
    config.innerRefineIterations = 30;
    config.outerVmag = 1:0.5:10;
    config.outerTheta = deg2rad(-180:10:180);
    config.outerPhi = 0;
    config.innerNu = linspace(0, 2*pi, 120);
    config.fixedLp = [];
end

function tofs = defaultGalileoTofs(planets, jd2k0, N, M)
    mu = GetBodyProps('Sun');
    body = planets{3};
    jd_earth1 = jd2k0 + 115 + 301;
    [sma_planet, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(body, jd_earth1);
    T_p = 2*pi*sqrt(sma_planet^3 / mu);
    tofs = [115, 301, T_p*M/86400, 1094];
    if N ~= 1
        tofs(3) = T_p*M/86400;
    end
end

function out = mergeDefaults(in, defaults)
    out = defaults;
    names = fieldnames(in);
    for idx = 1:numel(names)
        out.(names{idx}) = in.(names{idx});
    end
end

function planets = parsePlanetSequence(value)
    value = char(value);
    value = strrep(value, '->', ',');
    value = strrep(value, ';', ',');
    raw = strsplit(value, ',');
    planets = raw(~cellfun(@isempty, strtrim(raw)));
    for idx = 1:numel(planets)
        planets{idx} = strtrim(planets{idx});
    end
end

function values = normalizeResonanceVector(values)
    values = values(:)';
end

function apsisFlag = inferApsisFlag(N, M)
    if M > N
        apsisFlag = 1;
    elseif N > M
        apsisFlag = 0;
    else
        apsisFlag = 1;
    end
end

function resFlag = inferResFlag(resonantIndex, nTransfers)
    if resonantIndex == 1
        resFlag = 1;
    elseif resonantIndex == nTransfers
        resFlag = 3;
    else
        resFlag = 2;
    end
end

function config = computeAdjacentLambertContext(config)
    % Evaluate only non-resonant Lambert arcs. The selected resonant leg is
    % intentionally left open so ROGUI can inject a user-selected vinf_out.
    mr = 0;
    lp = 0;
    nTransfers = numel(config.planets) - 1;

    for idx = 1:nTransfers
        if config.is_vilm(idx)
            continue;
        end

        dnu = DeltaNu3(config.r(idx, :), config.r(idx + 1, :), 1);
        lw = double(dnu > pi);
        [config.vd(idx, :), config.va(idx, :)] = Lambert(config.r(idx, :), config.r(idx + 1, :), config.tofs(idx)*config.days2secs, config.mu_sun, lw, mr, lp);
    end

    % The resonant cost model is local: it needs the incoming excess
    % velocity from the previous leg and/or the required outgoing excess
    % velocity into the next leg, depending on where the resonance sits in
    % the sequence.
    idx = config.resonantIndex;
    if config.res_flag == 1
        config.vinfi = [0, 0, 0];
        config.vinfo = config.vd(idx + 1, :) - config.v(idx + 1, :);
    elseif config.res_flag == 3
        config.vinfi = config.va(idx - 1, :) - config.v(idx, :);
        config.vinfo = [0, 0, 0];
    else
        config.vinfi = config.va(idx - 1, :) - config.v(idx, :);
        config.vinfo = config.vd(idx + 1, :) - config.v(idx + 1, :);
    end
end
