%% Test script: compute and plot a simple MGA trajectory

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'AGA')));

% MGA Sequence Galileo (E-V-E-E-J), launch Oct 18 1989
% planets = {'Earth', 'Venus', 'Earth', 'Earth', 'Jupiter'};
% jd2k0 = -3727; % Oct 18, 1989

% JUICE
planets = {'Earth', 'Earth', 'Venus', 'Earth', 'Earth', 'Jupiter'};
jd2k0 = 8187;

% planets = {'Earth', 'Earth', 'Mars', 'Earth', 'Earth', 'Jupiter'};
% jd2k0 = 8493;



% N = [NaN,   NaN, 1,   NaN];
% M = [NaN,   NaN, 2,   NaN];

N = [1,   NaN, NaN, 1,   NaN];
M = [1,   NaN, NaN, 2,   NaN];

mu = GetBodyProps('Sun');
body = 'Earth';
% jd_earth1 = jd2k0 + 115 + 301;
jd_earth1 = jd2k0 + 175 + 332;
% jd_earth1 = jd2k0 + 823 + 571;
% [sma_planet, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(body, jd2k0);
% T_p0 = 2*pi*sqrt(sma_planet^3 / mu); % Earth osrbital period [s]
[sma_planet, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(body, jd_earth1);
T_p = 2*pi*sqrt(sma_planet^3 / mu); % Earth orbital period [s]

% tofs = [115, 301, T_p*M(3)/86400, 1094];
% tofs = [364, 823, 571, T_p*M(4)/86400, 1066]; % E->V, V->E1, E1->E2 (VILM), E2->J
tofs = [367, 175, 332, T_p*M(4)/86400, 1185];


% Get planet strings
seq = Planets2Seq(planets);
lPlanet = length(planets);

% Basic parameters
AU = 1.49597871E8; % Astronomical Unit (AU) [km]
n = 100; % Orbit resolution

%% Compute Lambert MGA
opts = struct('N_mag',8, 'N_theta',12, 'N_phi',1, 'N_points_nu',100, 'N_refine_nu',120,'algorithm','grid+fmin');

% Multi-Gravity Assist
[ jd2k, r, v, vd, va, rpga, dvga, dvdsm, vilm_arcs ] = MGA2_PGA2 ( planets, jd2k0, tofs, N, M, opts);

% Initial/Final Delta-V
dvd = norm(vd(1,:) - v(1,:)); % Departure DeltaV
dva = norm(va(end,:) - v(end,:)); % Arrival DeltaV

% GA altitude
[ ~, vmr ] = arrayfun(@GetBodyProps,planets(2:end-1));
hpga = rpga(:) ./ vmr(:);

%% Compute orbits

% Planet orbits
orbP = cell(lPlanet,1);
for k=1:lPlanet
    orbP{k} = GetBodyOrbit(planets{k},jd2k(k),mu,n,1);
end

% Spacecraft orbit
% Identify resonant legs from vilm_arcs
is_vilm = false(lPlanet-1, 1);
for k = 1:(lPlanet-1)
    is_vilm(k) = ~isempty(vilm_arcs(k).body);
end

orbS = cell(lPlanet-1,1);
for k=1:(lPlanet-1)
    if ~is_vilm(k)
        orbS{k} = ICF2Arc(r(k,:),vd(k,:),r(k+1,:),va(k,:),mu,100);
    end
end

orbit_res_arr = repmat(struct('rdsm',[NaN NaN NaN],'bdsm',[],'adsm',[], 'nu_deg',NaN,'revs',NaN,'dV_DSM',NaN), lPlanet-1, 1);
for k = 1:(lPlanet-1)
    if ~is_vilm(k)
        continue;
    end
    arcs = vilm_arcs(k);

    [~, nu_best, revs_best, ~, r_m_arc, v_m_minus_arc, v_m_plus_arc, ~, va_arc] = findOptimalDSMParameters(arcs.vinf_out, arcs.body, arcs.jd0, arcs.T_p, arcs.N, arcs.M, arcs.apsis_flag, arcs.mu_sun, arcs.search_nu, true);

    [r_p0, v_p0] = GetBodyICF(arcs.body, arcs.jd0, arcs.mu_sun, 0);
    sec2days = 1/86400;
    jd_f = arcs.jd0 + (arcs.M * arcs.T_p) * sec2days;
    [r_pf, ~] = GetBodyICF(arcs.body, jd_f, arcs.mu_sun, 0);

    v_sc0 = v_p0 + arcs.vinf_out;
    mr_lambert = arcs.N - revs_best - 1;

    [r_arc1, r_arc2, rdsm] = generateResonantTrajectoryPoints(r_p0, v_sc0, r_m_arc, v_m_minus_arc, v_m_plus_arc, r_pf, va_arc, revs_best, mr_lambert, arcs.mu_sun);

    orbit_res_arr(k) = struct('rdsm',rdsm,'bdsm',r_arc1,'adsm',r_arc2, 'nu_deg',rad2deg(nu_best),'revs',revs_best,'dV_DSM',dvdsm(k));
end

%% Plot MGA trajectory

% Reset figure
hold off;
plot3(0,0,0);

% Figure parameters
hold on; % Plot several data arrays
box on; % Box the canvas
axis equal; % Equal axis
view(0,90); % View camera
set(gca,'FontSize',12); % Font size

% Calendar date
[ year, month, day, ~, ~, ~ ] = JD2Cal(jd2k0,'J2000');

% Build r_pi string (filter NaN from VILM nodes without a flyby)
hpga_valid = hpga(~isnan(hpga));
if isempty(hpga_valid)
    rpi_str = 'r_\pi: N/A (VILM)';
else
    rpi_str = ['r_\pi: ',sprintf(' %.1f',hpga_valid),' R_{planet}'];
end

% Title
title({sprintf('MGA 2D PGA 2D | %s',seq);...
    [sprintf('T0: %04.0f-%02.0f-%02.0f',year,month,day),...
    ' | \Deltat:',sprintf(' %.0f',tofs),' days (',...
    sprintf('%.0f',sum(tofs)/365.25),'y)'];...
    ['\DeltaV:',sprintf(' %.0f',dvd*1E3),',',...
    sprintf(' %.0f,', dvga*1E3),...
    sprintf(' %.0f',dva*1E3),' m/s'];...
    ['\DeltaV_{DSM}:',sprintf(' %.0f',dvdsm(dvdsm>0)*1E3),' m/s'];...
    rpi_str});

% Plot planets
ph = cell(length(orbP),1);
col = cell(length(orbP),1); % Save colors
for k=1:length(orbP)
    is_first = ~any(strcmp(planets(1:k-1), planets{k}));
    if is_first
        trj = [orbP{k}(:,1),orbP{k}(:,2),orbP{k}(:,3)]'/AU;
        plot3(trj(1,:),trj(2,:),trj(3,:),'k--','LineWidth',0.5);
    end
    ph{k} = plot3(r(k,1)/AU,r(k,2)/AU,r(k,3)/AU,'o','LineWidth',1);
    col{k} = get(ph{k},'Color'); % Save assigned color
end

% Plot spacecraft
sh = plot3(0,0,0,'k-','LineWidth',0.5); % Plot spacecraft marker
for k=1:length(orbS) % Plot spacecraft orbit
    if is_vilm(k)
        % Plot two arcs + DSM point
        % Arc 1 (pre-DSM)
        trj1 = orbit_res_arr(k).bdsm' / AU;
        plot3(trj1(1,:), trj1(2,:), trj1(3,:), 'k-', 'LineWidth', 0.5);
        % Arc 2 (post-DSM)
        trj2 = orbit_res_arr(k).adsm' / AU;
        plot3(trj2(1,:), trj2(2,:), trj2(3,:), 'k-', 'LineWidth', 0.5);
        % DSM point
        dsm_h = plot3(orbit_res_arr(k).rdsm(1)/AU, orbit_res_arr(k).rdsm(2)/AU, orbit_res_arr(k).rdsm(3)/AU, ...
                       'rx','LineWidth',1,'MarkerSize',7);
    else
        trj = [orbS{k}(:,1),orbS{k}(:,2),orbS{k}(:,3)]'/AU; % Trajectory
        plot3(trj(1,:),trj(2,:),trj(3,:),'k-','LineWidth',0.5); % Plot Trj
    end
end

% Plot Sun marker
oh = plot3(0,0,0,'kx','LineWidth',1,'MarkerSize',1); % Sun

% Legend
seen = cell(length(orbP), 1);
ph_unique = cell(length(orbP), 1);
planets_unique = cell(length(orbP), 1);
n_unique = 0;

for k = 1:length(orbP)
    if ~any(strcmp(seen(1:n_unique), planets{k}))
        n_unique = n_unique + 1;
        seen{n_unique} = planets{k};
        ph_unique{n_unique} = ph{k};
        planets_unique{n_unique} = planets{k};
    else
        idx_first = find(strcmp(seen(1:n_unique), planets{k}), 1);
        set(ph{k}, 'Color', col{idx_first});
    end
end

ph_unique = ph_unique(1:n_unique);
planets_unique = planets_unique(1:n_unique);

legend([oh, ph_unique{:}, sh, dsm_h], 'Sun', planets_unique{:}, 'Spacecraft', 'DSM', 'Location', 'NorthEastOutside');


