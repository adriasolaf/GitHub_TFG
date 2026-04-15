function [jd2k, r, v, vd, va, rpga, dvga, dvdsm, orbit_res] = MGA2_PGA2(planets, jd2k0, tofs, N, M)
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
%   orbit_res: struct with resonant trajectory arc data
%
% References:
%   [-] n/a
%
% See also:
%   Lambert, GA_PGA2_Rp, OptimitzationVILM
%
% David de la Torre Sangra (original MGA framework)
% Adria Sola Foixench (VILM integration)
% January 2015 / April 2026

if nargin < 4
    N = 1; M = 1;
end

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
orbit_res = struct();

% Dates of encounters with planets
jd2k(1) = jd2k0;
for i = 2:lplanets
    jd2k(i) = jd2k(i-1) + tofs(i-1);
end

% Planet state vectors at each respective encounter
for i = 1:lplanets
    [r(i,:), v(i,:)] = GetBodyICF(planets{i}, jd2k(i), mu, 1);
end

% Detect resonant legs (same planet at consecutive encounters)
is_vilm = false(ltransfers, 1);
for i = 1:ltransfers
    if strcmp(planets{i}, planets{i+1})
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
vilm_count = 0;

for i = 1:ltransfers
    if ~is_vilm(i)
        continue;
    end
    vilm_count = vilm_count + 1;

    % Determine resonance position in the sequence
    if i == 1
        res_flag = 1; % Departure leg: dV_GA1 = 0
    elseif i == ltransfers
        res_flag = 3; % Arrival leg: dV_GA2 = 0
    else
        res_flag = 2; % Intermediate: full VILM cost
    end

    % Assign N and M
    if isscalar(N)
        n_val = N;
    else
        n_val = N(vilm_count);
    end

    if isscalar(M)
        m_val = M;
    else
        m_val = M(vilm_count);
    end

    % Automatic apsis flag based on resonance ratio
    if m_val > n_val
        apsis_flag = 1; % Outer orbit: DSM near apoapsis
    elseif n_val > m_val
        apsis_flag = 0; % Inner orbit: DSM near periapsis
    else
        apsis_flag = 1; % 1:1 resonance: default to apoapsis
    end

    % Extract v-infinity vectors from adjacent Lambert arcs
    if res_flag == 1
        vinfi = [0, 0, 0]; % No inbound arc
        vinfo = vd(i+1,:) - v(i+1,:);
    elseif res_flag == 3
        vinfi = va(i-1,:) - v(i,:);
        vinfo = [0, 0, 0]; % No outbound arc
    else
        vinfi = va(i-1,:) - v(i,:);
        vinfo = vd(i+1,:) - v(i+1,:);
    end

    p_name = planets{i};

    % Optimize the full VILM sequence
    [vout, ~, dV_GA1, dV_DSM, dV_GA2, ~, va_arr, ~, orbit_res, rp_GA1, rp_GA2] = OptimitzationVILM(p_name, jd2k(i), vinfi, vinfo, n_val, m_val, apsis_flag, mu, res_flag);

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
