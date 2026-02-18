%% Test script: compute and plot a simple MGA trajectory

% MGA Sequence
planets = {'Earth','Jupiter','Saturn'}; jd2k0 = -8169; tofs = [688,778]; % EJS

% Get planet strings
seq = Planets2Seq(planets);
lPlanet = length(planets);

% Basic parameters
mu = GetBodyProps('Sun'); % Standard gravitational parameter, Sun [km3/s2]
AU = 1.49597871E8; % Astronomical Unit (AU) [km]
n = 100; % Orbit resolution

%% Compute Lambert MGA

% Multi-Gravity Assist
[ jd2k, r, v, vd, va, rpga, dvga ] = MGA2_PGA2 ( planets, jd2k0, tofs );

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
orbS = cell(lPlanet-1,1);
for k=1:(lPlanet-1)
    orbS{k} = ICF2Arc(r(k,:),vd(k,:),r(k+1,:),va(k,:),mu,100);
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

% Title
title({sprintf('MGA 2D PGA 2D | %s',seq);...
    [sprintf('T0: %04.0f-%02.0f-%02.0f',year,month,day),...
    ' | \Deltat:',sprintf(' %.0f',tofs),' days (',...
    sprintf('%.0f',sum(tofs)/365.25),'y)'];...
    ['\DeltaV:',sprintf(' %.0f',dvd*1E3),', ',...
    sprintf(' %.0f',dvga*1E3),', ',...
    sprintf(' %.0f',dva*1E3),' m/s'];...
    ['r_\pi: ',sprintf(' %.0f',hpga),' R_{planet}']});

% Plot planets
ph = cell(length(orbP),1);
for k=1:length(orbP)

    % Plot planet(i) orbit
    ce = strfind(seq,seq(k)); % Index of current encounter inside sequence
    if ce(1)==k % Only one orbit for each planet, regardless of encounters
        trj = [orbP{k}(:,1),orbP{k}(:,2),orbP{k}(:,3)]'/AU; % Trajectory
        plot3(trj(1,:),trj(2,:),trj(3,:),'k--','LineWidth',0.5); % Plot Trj
    end

    % Plot encounter marker
    ph{k} = plot3(r(k,1)/AU,r(k,2)/AU,r(k,3)/AU,'o','LineWidth',1);

end

% Plot spacecraft
sh = plot3(0,0,0,'k-','LineWidth',0.5); % Plot spacecraft marker
for k=1:length(orbS) % Plot spacecraft orbit
    trj = [orbS{k}(:,1),orbS{k}(:,2),orbS{k}(:,3)]'/AU; % Trajectory
    plot3(trj(1,:),trj(2,:),trj(3,:),'k-','LineWidth',0.5); % Plot Trj
end

% Plot Sun marker
oh = plot3(0,0,0,'kx','LineWidth',1,'MarkerSize',1); % Sun

% Legend
legend([oh,ph{:},sh],'Sun',planets{:},'Spacecraft',...
    'Location','NorthEastOutside');
