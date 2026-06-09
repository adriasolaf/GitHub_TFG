function h = plotOuterSearchMap(ax, outerMap, selectedPoint, metricName, refinedPoint)
%plotOuterSearchMap Plot outer search map
%   Plots the outer ROGUI search table over theta, phi, and vmag. Angles in
%   the table are radians and are displayed in degrees.
%
% Inputs:
%   ax: target MATLAB axes
%   outerMap: outer search result table
%   selectedPoint: selected outer row or empty
%   metricName: table variable used as color metric
%
% Outputs:
%   h: plot graphics handle structure
%
% Example:
%   [ h ] = plotOuterSearchMap ( ax, outerMap, bestRow, 'totalDv' );
%
% References:
%   [-]
%
%May 2026

    if nargin < 1 || ~isgraphics(ax, 'axes')
        ax = gca;
    end
    if nargin < 2
        outerMap = table();
    end
    if nargin < 3
        selectedPoint = [];
    end
    if nargin < 4 || isempty(metricName)
        metricName = 'totalDv';
    end
    if nargin < 5
        refinedPoint = [];
    end

    cla(ax);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');
    h = struct('map', gobjects(0), 'selected', gobjects(0), 'refined', gobjects(0), 'colorbar', gobjects(0));

    if isempty(outerMap) || ~istable(outerMap) || ~all(ismember({'theta','phi','vmag'}, outerMap.Properties.VariableNames))
        showEmptyAxes(ax, 'Outer search map unavailable');
        hold(ax, 'off');
        return;
    end
    if ~any(strcmp(metricName, outerMap.Properties.VariableNames))
        metricName = 'totalDv';
    end

    values = outerMap.(metricName);
    if any(strcmp('isFeasible', outerMap.Properties.VariableNames))
        values(~outerMap.isFeasible) = NaN;
    end

    thetaDeg = rad2deg(outerMap.theta);
    phiDeg = rad2deg(outerMap.phi);
    phiValues = unique(round(phiDeg * 1e10) / 1e10);
    is3d = numel(phiValues) > 1;

    if is3d
        % A true 3D search cannot be represented as a single surface without
        % choosing a slice, so use branch-colored samples.
        valid = isfinite(thetaDeg) & isfinite(phiDeg) & isfinite(outerMap.vmag) & isfinite(values);
        if any(valid)
            h.map = scatter3(ax, thetaDeg(valid), phiDeg(valid), outerMap.vmag(valid), 30, values(valid), 'filled');
            xlabel(ax, '\theta [deg]');
            ylabel(ax, '\phi [deg]');
            zlabel(ax, '|v_\infty| [km/s]');
            view(ax, 3);
        end
    else
        % The common case is a fixed phi slice. Convert the sampled table
        % back into a dense grid so contourf can show topology more clearly
        % than a point cloud.
        [thetaGrid, vmagGrid, valueGrid] = buildGrid(thetaDeg, outerMap.vmag, values);
        if ~isempty(valueGrid) && any(isfinite(valueGrid(:)))
            [~, h.map] = contourf(ax, thetaGrid, vmagGrid, valueGrid, 20, 'LineStyle', 'none');
        end
        xlabel(ax, '\theta [deg]');
        ylabel(ax, '|v_\infty| [km/s]');
    end

    if isempty(h.map) || ~isgraphics(h.map(1))
        showEmptyAxes(ax, 'Outer search map has no finite values');
        hold(ax, 'off');
        return;
    end

    colormap(ax, 'parula');
    h.colorbar = colorbar(ax);
    ylabel(h.colorbar, metricName, 'Interpreter', 'none');
    title(ax, 'Outer search map');

    marker = selectedMarker(selectedPoint, outerMap);
    if ~isempty(marker)
        if is3d
            h.selected = plot3(ax, marker(1), marker(2), marker(3), 'rx', 'LineWidth', 2, 'MarkerSize', 11);
        else
            h.selected = plot(ax, marker(1), marker(3), 'rx', 'LineWidth', 2, 'MarkerSize', 11);
        end
    end

    
    refMarker = refinedMarker(refinedPoint);
    if ~isempty(refMarker)
        if is3d
            h.refined = plot3(ax, refMarker(1), refMarker(2), refMarker(3), 'p', 'MarkerFaceColor', [0 0.7 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 14, 'LineWidth', 1);
        else
            h.refined = plot(ax, refMarker(1), refMarker(3), 'p', 'MarkerFaceColor', [0 0.7 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 14, 'LineWidth', 1);
        end
    end

    hold(ax, 'off');
end

function marker = refinedMarker(refinedPoint)
    marker = [];
    if isempty(refinedPoint)
        return;
    end
    if isstruct(refinedPoint) && isfield(refinedPoint, 'theta')
        marker = [rad2deg(refinedPoint.theta), rad2deg(refinedPoint.phi), refinedPoint.vmag];
    end
end

function [X, Y, Z] = buildGrid(x, y, values)
    % Reconstruct a rectangular grid from table rows. Missing or infeasible
    % samples remain NaN so contourf leaves holes instead of inventing data.
    X = [];
    Y = [];
    Z = [];
    valid = isfinite(x) & isfinite(y) & isfinite(values);
    if ~any(valid)
        return;
    end

    xUnique = unique(x(valid));
    yUnique = unique(y(valid));
    if numel(xUnique) < 2 || numel(yUnique) < 2
        return;
    end

    Z = NaN(numel(yUnique), numel(xUnique));
    for idx = find(valid(:))'
        ix = find(xUnique == x(idx), 1);
        iy = find(yUnique == y(idx), 1);
        Z(iy, ix) = values(idx);
    end
    [X, Y] = meshgrid(xUnique, yUnique);
end

function marker = selectedMarker(selectedPoint, outerMap)
    % If no explicit selection is available, mark the best feasible sample
    % so users can immediately see the optimizer's current preference.
    marker = [];
    if isempty(selectedPoint)
        feasible = isfinite(outerMap.totalDv);
        if any(strcmp('isFeasible', outerMap.Properties.VariableNames))
            feasible = feasible & outerMap.isFeasible;
        end
        if ~any(feasible)
            return;
        end
        rows = outerMap(feasible, :);
        [~, idx] = min(rows.totalDv);
        selectedPoint = rows(idx, :);
    end

    if istable(selectedPoint)
        marker = [rad2deg(selectedPoint.theta(1)), rad2deg(selectedPoint.phi(1)), selectedPoint.vmag(1)];
    elseif isstruct(selectedPoint) && isfield(selectedPoint, 'theta')
        marker = [rad2deg(selectedPoint.theta), rad2deg(selectedPoint.phi), selectedPoint.vmag];
    elseif isnumeric(selectedPoint) && numel(selectedPoint) >= 3
        marker = selectedPoint(1:3);
    end
end

function showEmptyAxes(ax, message)
    text(ax, 0.5, 0.5, message, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'Color', [0.35 0.35 0.35]);
    xlabel(ax, '\theta [deg]');
    ylabel(ax, '|v_\infty| [km/s]');
    title(ax, 'Outer search map');
end
