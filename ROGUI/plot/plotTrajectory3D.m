function h = plotTrajectory3D(ax, r, orbS, orbitRes, planets, orbP)
% plotTrajectory3D Plot planet encounters, spacecraft arcs, and DSM points.
%
% Usage:
%   h = plotTrajectory3D(ax, r, orbS, orbitRes, planets, orbP)
%   h = plotTrajectory3D(ax, dataStruct)
%
% Position inputs are expected in km and are plotted in AU. The plotted
% z-coordinate is visually exaggerated by zScale for inspection only; this
% does not modify the trajectory data returned by evaluation helpers.

if nargin < 1 || ~isgraphics(ax, 'axes')
    ax = gca;
end
if nargin < 2
    r = [];
end
if nargin < 3
    orbS = {};
end
if nargin < 4
    orbitRes = [];
end
if nargin < 5
    planets = {};
end
if nargin < 6
    orbP = {};
end

if isstruct(r)
    data = r;
    if isfield(data, 'plot') && isstruct(data.plot)
        orbP = getFirstField(data.plot, {'planetOrbits'}, {});
        orbS = getFirstField(data.plot, {'transferArcs'}, {});
        orbitRes = getFirstField(data.plot, {'resonant','resonantArcs'}, []);
    else
        orbP = getFirstField(data, {'orbP','planetOrbits'}, {});
        orbS = getFirstField(data, {'orbS','spacecraftArcs','trajectoryArcs','transferArcs'}, {});
        orbitRes = getFirstField(data, {'orbitRes','orbit_res_arr','resonantArcs','resonant'}, []);
    end
    planets = getFirstField(data, {'planets','bodies'}, {});
    r = getFirstField(data, {'r','encounters','planetPositions'}, []);
end

% Preserve the user's current camera between selection updates. A blank
% axes defaults to the mission-design top view requested for ROGUI.
if isempty(ax.Children)
    oldView = [0 90];
else
    [az, el] = view(ax);
    oldView = [az, el];
end

cla(ax);
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
view(ax, oldView);

AU = 1.49597871E8;
zScale = 10; % Visual-only vertical exaggeration.
h = struct('sun', gobjects(0), 'planetOrbits', gobjects(0), ...
    'planets', gobjects(0), 'spacecraft', gobjects(0), 'dsm', gobjects(0));

try
    h.sun = plot3(ax, 0, 0, 0, 'kx', 'LineWidth', 1, 'MarkerSize', 6);
    h.planetOrbits = plotPlanetOrbits(ax, orbP, AU, zScale);
    h.planets = plotEncounters(ax, r, planets, AU, zScale);
    h.spacecraft = plotSpacecraftArcs(ax, orbS, orbitRes, AU, zScale);
    h.dsm = plotDsmPoints(ax, orbitRes, AU, zScale);

    xlabel(ax, 'x [AU]');
    ylabel(ax, 'y [AU]');
    zlabel(ax, 'z [AU] x10');
    title(ax, 'Trajectory');
    daspect(ax, [1 1 1]);
    view(ax, oldView);
catch
    cla(ax);
    text(ax, 0.5, 0.5, 'Trajectory could not be plotted', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'Color', [0.35 0.35 0.35]);
    view(ax, oldView);
end

if isempty([h.planetOrbits(:); h.planets(:); h.spacecraft(:); h.dsm(:)])
    text(ax, 0.5, 0.5, 'Trajectory unavailable', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'Color', [0.35 0.35 0.35]);
end

hold(ax, 'off');
end

function handles = plotPlanetOrbits(ax, orbP, AU, zScale)
handles = gobjects(0);
if ~iscell(orbP)
    orbP = {orbP};
end
for k = 1:numel(orbP)
    pts = normalizePoints(orbP{k});
    if isempty(pts)
        continue;
    end
    handles(end+1,1) = plot3(ax, pts(:,1)/AU, pts(:,2)/AU, zScale*pts(:,3)/AU, ...
        'k--', 'LineWidth', 0.5); %#ok<AGROW>
end
end

function handles = plotEncounters(ax, r, planets, AU, zScale)
handles = gobjects(0);
pts = normalizePoints(r);
if isempty(pts)
    return;
end

for k = 1:size(pts, 1)
    handles(end+1,1) = plot3(ax, pts(k,1)/AU, pts(k,2)/AU, zScale*pts(k,3)/AU, ...
        'o', 'LineWidth', 1, 'MarkerSize', 5); %#ok<AGROW>
    if iscell(planets) && k <= numel(planets) && ~isempty(planets{k})
        text(ax, pts(k,1)/AU, pts(k,2)/AU, zScale*pts(k,3)/AU, [' ', planets{k}], ...
            'FontSize', 8);
    end
end
end

function handles = plotSpacecraftArcs(ax, orbS, orbitRes, AU, zScale)
handles = gobjects(0);
if ~iscell(orbS)
    orbS = {orbS};
end
for k = 1:numel(orbS)
    pts = normalizePoints(orbS{k});
    if ~isempty(pts)
        handles(end+1,1) = plot3(ax, pts(:,1)/AU, pts(:,2)/AU, zScale*pts(:,3)/AU, ...
            'k-', 'LineWidth', 0.8); %#ok<AGROW>
    end
end

if isempty(orbitRes) || ~isstruct(orbitRes)
    return;
end
for k = 1:numel(orbitRes)
    if isfield(orbitRes(k), 'bdsm')
        pts = normalizePoints(orbitRes(k).bdsm);
        if ~isempty(pts)
            handles(end+1,1) = plot3(ax, pts(:,1)/AU, pts(:,2)/AU, zScale*pts(:,3)/AU, ...
                'k-', 'LineWidth', 0.8); %#ok<AGROW>
        end
    end
    if isfield(orbitRes(k), 'adsm')
        pts = normalizePoints(orbitRes(k).adsm);
        if ~isempty(pts)
            handles(end+1,1) = plot3(ax, pts(:,1)/AU, pts(:,2)/AU, zScale*pts(:,3)/AU, ...
                'k-', 'LineWidth', 0.8); %#ok<AGROW>
        end
    end
    if isfield(orbitRes(k), 'arc1')
        pts = normalizePoints(orbitRes(k).arc1);
        if ~isempty(pts)
            handles(end+1,1) = plot3(ax, pts(:,1)/AU, pts(:,2)/AU, zScale*pts(:,3)/AU, ...
                'k-', 'LineWidth', 0.8); %#ok<AGROW>
        end
    end
    if isfield(orbitRes(k), 'arc2')
        pts = normalizePoints(orbitRes(k).arc2);
        if ~isempty(pts)
            handles(end+1,1) = plot3(ax, pts(:,1)/AU, pts(:,2)/AU, zScale*pts(:,3)/AU, ...
                'k-', 'LineWidth', 0.8); %#ok<AGROW>
        end
    end
end
end

function handles = plotDsmPoints(ax, orbitRes, AU, zScale)
handles = gobjects(0);
if isempty(orbitRes) || ~isstruct(orbitRes)
    return;
end
for k = 1:numel(orbitRes)
    if ~isfield(orbitRes(k), 'rdsm')
        if isfield(orbitRes(k), 'rDsm')
            dsmPoint = orbitRes(k).rDsm;
        else
            continue;
        end
    else
        dsmPoint = orbitRes(k).rdsm;
    end
    pts = normalizePoints(dsmPoint);
    if isempty(pts)
        continue;
    end
    handles(end+1,1) = plot3(ax, pts(1,1)/AU, pts(1,2)/AU, zScale*pts(1,3)/AU, ...
        'rx', 'LineWidth', 1.2, 'MarkerSize', 7); %#ok<AGROW>
end
end

function pts = normalizePoints(pts)
if isempty(pts) || ~isnumeric(pts)
    pts = [];
    return;
end
if size(pts, 2) ~= 3 && size(pts, 1) == 3
    pts = pts';
end
if size(pts, 2) ~= 3
    pts = [];
    return;
end
pts = pts(all(isfinite(pts), 2), :);
end

function value = getFirstField(s, names, defaultValue)
value = defaultValue;
for k = 1:numel(names)
    if isfield(s, names{k}) && ~isempty(s.(names{k}))
        value = s.(names{k});
        return;
    end
end
end
