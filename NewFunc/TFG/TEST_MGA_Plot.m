%% Test script: compute and plot a simple MGA trajectory

% MGA Sequence
planets = {'Earth', 'Jupiter', 'Jupiter', 'Saturn'}; jd2k0 = -8169;

N=2; M=3;

tofs = [688,11.86*365.25636*M,5*778]; % EJS

% Get planet strings
seq = Planets2Seq(planets);
lPlanet = length(planets);

% Basic parameters
mu = GetBodyProps('Sun'); % Standard gravitational parameter, Sun [km3/s2]
AU = 1.49597871E8; % Astronomical Unit (AU) [km]
n = 100; % Orbit resolution

%% Compute Lambert MGA

% Multi-Gravity Assist
[ jd2k, r, v, vd, va, rpga, dvga, dvdsm, orbita_res ] = MGA2_PGA2 ( planets, jd2k0, tofs, N, M );

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
% Identify resonant legs (same planet at both ends)
is_vilm = false(lPlanet-1, 1);
for k = 1:(lPlanet-1)
    if strcmp(planets{k}, planets{k+1})
        is_vilm(k) = true;
    end
end

orbS = cell(lPlanet-1,1);
for k=1:(lPlanet-1)
    if ~is_vilm(k)
        orbS{k} = ICF2Arc(r(k,:),vd(k,:),r(k+1,:),va(k,:),mu,100);
    end
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
        trj1 = orbita_res.bdsm' / AU;
        plot3(trj1(1,:), trj1(2,:), trj1(3,:), 'k-', 'LineWidth', 0.5);
        % Arc 2 (post-DSM)
        trj2 = orbita_res.adsm' / AU;
        plot3(trj2(1,:), trj2(2,:), trj2(3,:), 'k-', 'LineWidth', 0.5);
        % DSM point
        dsm_h = plot3(orbita_res.rdsm(1)/AU, orbita_res.rdsm(2)/AU, orbita_res.rdsm(3)/AU, ...
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