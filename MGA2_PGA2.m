function [jd2k, r, v, vd, va, rpga, dvga, dvdsm, vilm_arcs] = MGA2_PGA2(planets, jd2k0, tofs, N, M, opts)
%   Performs a multi-gravity assist trajectory by means of:
%   - Powered Gravity Assist, 2D
%   - Lambert arcs
%   - Resonant VILM orbits
%   The algorithm relies on the patched conics method.
%
% Inputs:
%   planets: encounter planets sequence (cell array of strings)
%   jd2k0: sequence initial date [days from J2000]
%   tofs: times of flight between sequence objects [days]
%   N: number of spacecraft revolutions (scalar or vector per VILM)
%   M: number of body revolutions (scalar or vector per VILM)
%   opts: struct with optimizer settings forwarded to optimizeOutgoingVInfinityVILM. Same opts is
%         applied to every VILM leg in the sequence.
%
% Outputs:
%   jd2k: time at each encounter [days from J2000]
%   r: heliocentric position of planet at each encounter [km]
%   v: heliocentric velocity of planet at each encounter [km/s]
%   vd: heliocentric departure velocity of s/c at each transfer [km/s]
%   va: heliocentric arrival velocity of s/c at each transfer [km/s]
%   dvga: Gravity assist DeltaV at each encounter [km/s]
%   rpga: Gravity assist periapsis radius at each encounter [km]
%   dvdsm: DSM DeltaV for resonant legs [km/s]
%   vilm_arcs: struct array with the inputs needed to rebuild each
%              resonant leg externally (one entry per transfer; non-VILM
%              entries are empty)
%
% References:
%   [-] n/a
%
% See also:
%   Lambert, GA_PGA2_Rp, optimizeOutgoingVInfinityVILM
%
% David de la Torre Sangra (original MGA framework)
% Adria Sola Foixench (Resonant integration)
% January 2015 / May 2026

% Constants
mu = GetBodyProps('Sun');
days2secs = 86400;

% Lambert configuration
mr = 0; % No multi-revolutions
lp = 0; % Short-period solutions

% Auxiliary magnitudes
lplanets = length(planets);
ltransfers = lplanets - 1;

% Preallocate arrays
jd2k = zeros(lplanets,1);
r = zeros(lplanets,3);
v = zeros(lplanets,3);
vd = zeros(ltransfers,3);
va = zeros(ltransfers,3);
dvga = zeros(ltransfers-1,1);
rpga = zeros(ltransfers-1,1);
dvdsm = zeros(ltransfers,1);
vilm_arcs = repmat(struct('vinf_out',[],'body','','jd0',[],'T_p',[], 'N',[],'M',[],'apsis_flag',[],'mu_sun',[],'search_nu',[], 'vinfi_req',[],'vinfo_req',[],'res_flag',[]), ltransfers, 1);

% Dates of encounters with planets
jd2k(1) = jd2k0;
for i = 2:lplanets
    jd2k(i) = jd2k(i-1) + tofs(i-1);
end

% Planet state vectors at each respective encounter
for i = 1:lplanets
    [r(i,:), v(i,:)] = GetBodyICF(planets{i}, jd2k(i), mu, 1);
end

% Build N and M arrays per length transfers.
[N_per_leg, flagN] = scalar2LengthN(N, ltransfers);
[M_per_leg, flagM] = scalar2LengthN(M, ltransfers);

if flagN ~= 0 || flagM ~= 0
    return;
end

% Detect resonant legs: same planet at consecutive encounters and
% tofs(i) matches the N:M resonant period (M * T_p) within tolerance.
is_vilm = false(ltransfers, 1);
for i = 1:ltransfers
    if ~strcmp(planets{i}, planets{i+1})
        continue;
    end
    m_val = M_per_leg(i);
    if isnan(m_val)
        continue;
    end

    % Orbital period of the resonant body
    [sma_p, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(planets{i}, jd2k(i));
    T_p_s = 2*pi * sqrt(sma_p^3 / mu);

    % Resonant tof
    tof_res = m_val * T_p_s / days2secs;

    tol_days = 2;
    if abs(tofs(i) - tof_res) < tol_days
        is_vilm(i) = true;
    end
end


% Lambert transfer arcs (non-resonant legs)
for i = 1:ltransfers
    if is_vilm(i)
        continue; % Resonant leg, handled separately
    end

    % Transfer angle for prograde motion
    dnu = DeltaNu3(r(i,:), r(i+1,:), 1);

    if dnu > pi
        lw = 1;
    else
        lw = 0;
    end

    % Lambert arc from planet i to planet i+1
    [vd(i,:), va(i,:)] = Lambert(r(i,:), r(i+1,:), tofs(i)*days2secs, mu, lw, mr, lp);
end


% Resonant VILM orbits
for i = 1:ltransfers
    if ~is_vilm(i)
        continue;
    end

    % Determine resonance position in the sequence
    if i == 1
        res_flag = 1; % Departure leg: dV_GA1 = 0
    elseif i == ltransfers
        res_flag = 3; % Arrival leg: dV_GA2 = 0
    else
        res_flag = 2; % Intermediate: full VILM cost
    end

    % N and M for this leg (indexed by transfer index, not by counter)
    n_val = N_per_leg(i);
    m_val = M_per_leg(i);

    % Automatic apsis flag based on resonance ratio
    if m_val > n_val
        apsis_flag = 1; % Outer transfer: anchor the DSM scan at apoapsis
    elseif n_val > m_val
        apsis_flag = 0; % Inner transfer: anchor the DSM scan at periapsis
    else
        apsis_flag = 1; % 1:1 resonance: default scan anchor at apoapsis
    end

    % Extract v-infinity vectors from adjacent Lambert arcs
    if res_flag == 1
        vinfi_req = [0, 0, 0]; % No inbound arc
        vinfo_req = vd(i+1,:) - v(i+1,:);
    elseif res_flag == 3
        vinfi_req = va(i-1,:) - v(i,:);
        vinfo_req = [0, 0, 0]; % No outbound arc
    else
        vinfi_req = va(i-1,:) - v(i,:);
        vinfo_req = vd(i+1,:) - v(i+1,:);
    end

    p_name = planets{i};

    % Optimize the full VILM sequence
    [vout, ~, dV_GA1, dV_DSM, dV_GA2, ~, va_arr, ~, vilm_arc_i, rp_GA1, rp_GA2] = optimizeOutgoingVInfinityVILM(p_name, jd2k(i), vinfi_req, vinfo_req, n_val, m_val, apsis_flag, mu, res_flag, opts);

    vilm_arcs(i) = vilm_arc_i;

    % Store results
    vd(i,:) = v(i,:) + vout;
    va(i,:) = va_arr;
    dvdsm(i) = dV_DSM;

    % Assign gravity assist costs to the corresponding encounter index
    if res_flag == 1
        dvga(i) = dV_GA2; rpga(i) = rp_GA2;
    elseif res_flag == 2
        dvga(i-1) = dV_GA1; rpga(i-1) = rp_GA1;
        dvga(i) = dV_GA2; rpga(i) = rp_GA2;
    elseif res_flag == 3
        dvga(i-1) = dV_GA1; rpga(i-1) = rp_GA1;
    end
end


% Standard powered gravity assists (non-VILM encounters)
for i = 1:(ltransfers-1)
    if is_vilm(i) || is_vilm(i+1)
        continue; % GA already computed by the VILM optimizer
    end

    vinfi_ga = va(i,:) - v(i+1,:);
    vinfo_ga = vd(i+1,:) - v(i+1,:);
    vinfin_ga = norm(vinfi_ga);
    vinfon_ga = norm(vinfo_ga);

    % Deflection angle between inbound and outbound v-infinity
    cos_delta = dot(vinfi_ga, vinfo_ga) / (vinfin_ga * vinfon_ga);
    delta = acos(cos_delta);

    [mu_planet, vmr] = GetBodyProps(planets{i+1});
    [dvga(i), rpga(i)] = GA_PGA2_Rp(vinfin_ga, vinfon_ga, delta, mu_planet);

    vmr_safety = 1.05 * vmr;
    if rpga(i) < vmr_safety
        rpga(i) = vmr_safety;
        [dvga(i)] = GA_PGA2_Vinfo(vinfin_ga, vinfon_ga, delta, vmr_safety, mu_planet);
    end
end

end



