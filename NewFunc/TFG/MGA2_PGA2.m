function [jd2k, r, v, vd, va, rpga, dvga, dvdsm, orbita_res] = MGA2_PGA2(planets, jd2k0, tofs, N, M)
%   Performs a multi-gravity assist trajectory by means of:
%   - Powered Gravity Assist, 2D
%   - Lambert arcs
%   The algorithm relies on the patched conics method.
%
% Inputs:
%   planets: encounter planets sequence
%   jd2k0: sequence initial date [days from J2000]
%   tofs: times of flight between sequence objects [days]
%
% Outputs:
%   jd2k: time at each encounter [days from J2000]
%   r: heliocentric position of planet at each encounter [km]
%   v: heliocentric velocity of planet at each encounter [km/s]
%   vd: heliocentric departure velocity of s/c at each transfer [km/s]
%   va: heliocentric arrival velocity of s/c at each transfer [km/s]
%   dvga: Gravity assist DeltaV manoeure at each encounter [km/s]
%   rpga: Gravity assist periapsis radius at each encounter [km]
%
% Note: indices are considered as follows:
%   |
%   |-- Sequence: EMJ (Earth->Mars->Jupiter)
%   |
%   |-> Encounter i=1 (departure from Earth)
%   |   |-> jd2k(1), r(1,:), v(1,:)
%   |
%   |-> Transfer i=1 (transfer arc from Earth to Mars)
%   |   |-> vd(1,:), va(1,:)
%   |
%   |-> Encounter i=2 (GA with Mars)
%   |   |-> jd2k(2), r(2,:), v(2,:), dvga(1)
%   |
%   |-> Transfer i=2 (transfer arc from Mars to Jupiter)
%   |   |-> vd(2,:), va(2,:)
%   |
%   |-> Encounter i=3 (arrival to Jupiter)
%       |-> jd2k(3), r(3,:), v(3,:)
%
% Example:
%   <https://solarsystem.nasa.gov/multimedia/gallery/Voyager_Path.jpg>
%   seq = 'EJSUN'; t0 = -8169; tof = [688,778,1613,1309]; % Voyager 1
%   [ t, r, v, vd, va ] = MGA_PGA2 ( seq, t0, tof ); % Compute MGA
%   Voyager-1 c3: 105.5 (10.27km/s), Voyager-2 c3: 102.4 (10.12km/s)
%
% References:
%	[-] n/a
%
% See also:
%	Lambert, GA_PGA2
%
% David de la Torre Sangra
% January 2015

if nargin < 4
    N = 1; M = 1; 
end

% Constants
mu = GetBodyProps('Sun'); % Standard gravitational parameter (Sun) [km3/s2]
days2secs = 86400; % Days to seconds

% Lambert configuration
mr = 0; % No multi-revolutions
lp = 0; % Short-period solutions

% Auxiliary magnitudes
lplanets = length(planets); % Number of planets
ltransfers = lplanets - 1; % Number of transfers


% Preallocate arrays
jd2k = zeros(lplanets,1);
r = zeros(lplanets,3);
v = zeros(lplanets,3);
vd = zeros(ltransfers,3);
va = zeros(ltransfers,3);
dvga = zeros(ltransfers-1,1);
rpga = zeros(ltransfers-1,1);
dvdsm = zeros(ltransfers,1);
orbita_res = struct();

% Dates of encounters with planets
jd2k(1) = jd2k0; % Departure date
for i=2:lplanets
    jd2k(i) = jd2k(i-1) + tofs(i-1); % Encounter i
end

% Planet state vectors at each respective encounter
for i=1:lplanets
    [r(i,:), v(i,:)] = GetBodyICF(planets{i}, jd2k(i), mu, 1);
end

% Logic vector to find resonant orbits
is_vilm = false(ltransfers, 1);
for i = 1:ltransfers
    if strcmp(planets{i}, planets{i+1})
        is_vilm(i) = true;
    end
end

% Sequence of Lambert transfer arcs
for i=1:ltransfers % Iterate sequence of transfers

    if is_vilm(i)
        continue; % Not an pure Lambert arc
    end

    % Compute transfer angle, accounting for prograde motion
    dnu = DeltaNu3(r(i,:), r(i+1,:), 1); % Transfer angle

    if dnu > pi
        lw = 1; 
    else 
        lw = 0; 
    end % Long/short way tranfer

    % Compute Lambert arc from the current planet up to the next planet
    [vd(i,:), va(i,:)] = Lambert(r(i,:), r(i+1,:), tofs(i)*days2secs, mu, lw, mr, lp);
end

% Sequence of resonant orbits
vilm_count = 0; % Initialize counter in case of multiple VILM sequences

for i = 1:ltransfers
    if ~is_vilm(i)
        continue; % Not a resonant orbit
    end
    vilm_count = vilm_count + 1;
    
    % Automated auxiliar flags
    if i == 1
        res_flag = 1; 
    elseif i == ltransfers
        res_flag = 3; 
    else
        res_flag = 2; 
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
    
    % Automated apsis flags
    if m_val > n_val
        apsis_flag = 1; % Outer orbit -> DSM at Apoapsis
    elseif n_val > m_val
        apsis_flag = 0; % Inner orbit -> DSM at Periapsis
    else
        apsis_flag = 1; % 1:1 N:M
    end
    
    % Extract vinfi and vinfo from previous and next Lambert arcs
    if res_flag == 1
        vinfi = [0, 0, 0];
        vinfo = vd(i+1,:) - v(i+1,:);
    elseif res_flag == 3
        vinfi = va(i-1,:) - v(i,:);
        vinfo = [0, 0, 0];
    else
        vinfi = va(i-1,:) - v(i,:);
        vinfo = vd(i+1,:) - v(i+1,:);
    end
    
    p_name = planets{i};

    [vout, ~, dV_GA1, dV_DSM, dV_GA2, ~, va_arr, ~, orbita_res, rp_GA1, rp_GA2] = OptimitzationVILM(p_name, jd2k(i), vinfi, vinfo, n_val, m_val, apsis_flag, mu, res_flag);
        
    % Save results
    vd(i,:) = v(i,:) + vout;
    va(i,:) = va_arr;
    dvdsm(i) = dV_DSM;
    
    % Assign gravity assist impulses to the corresponding index
    if res_flag == 1
        dvga(i) = dV_GA2; rpga(i) = rp_GA2;
    elseif res_flag == 2
        dvga(i-1) = dV_GA1; rpga(i-1) = rp_GA1;
        dvga(i) = dV_GA2; rpga(i)   = rp_GA2;
    elseif res_flag == 3
        dvga(i-1) = dV_GA1; rpga(i-1) = rp_GA1;
    end
end


for i=1:(ltransfers-1) 
    if is_vilm(i) || is_vilm(i+1) % The entry or the departure transfer is already computed
        continue; 
    end

    vinfi = va(i,:) - v(i+1,:); 
    vinfo = vd(i+1,:) - v(i+1,:); 
    vinfin = norm(vinfi); 
    vinfon = norm(vinfo); 

    cos = dot(vinfi, vinfo) / (vinfin*vinfo);

    delta = acos(cos);

    [mu_planet, vmr] = GetBodyProps(planets{i+1});
    [dvga(i), rpga(i)] = GA_PGA2_Rp(vinfin, vinfon, delta, mu_planet);

    vmr_safety = 1.05 * vmr; 
    if rpga(i) < vmr_safety 
        rpga(i) = vmr_safety; 
        [dvga(i)] = GA_PGA2_Vinfo(vinfin, vinfon, delta, vmr_safety, mu_planet);
    end
end
end