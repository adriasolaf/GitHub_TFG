function h = plotInnerSearchMap(ax, innerMap, selectedPoint, metricName)
%plotInnerSearchMap Plot inner DSM search map
%   Plots branch-resolved DSM anomaly samples for one selected outgoing
%   v-infinity vector. nu is displayed in degrees.
%
% Inputs:
%   ax: target MATLAB axes
%   innerMap: inner search result table
%   selectedPoint: selected inner row or empty
%   metricName: table variable used on the vertical axis
%
% Outputs:
%   h: plot graphics handle structure
%
% Example:
%   [ h ] = plotInnerSearchMap ( ax, innerMap, bestRow, 'totalDv' );
%
% References:
%   [-]
%
%May 2026

    if nargin < 1 || ~isgraphics(ax, 'axes')
        ax = gca;
    end
    if nargin < 2
        innerMap = table();
    end
    if nargin < 3
        selectedPoint = [];
    end
    if nargin < 4 || isempty(metricName)
        metricName = 'totalDv';
    end

    cla(ax);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');
    h = struct('points', gobjects(0), 'selected', gobjects(0), 'colorbar', gobjects(0));

    if isempty(innerMap) || ~istable(innerMap) || ~any(strcmp('nu', innerMap.Properties.VariableNames))
        showEmptyAxes(ax, 'Inner search map unavailable');
        hold(ax, 'off');
        return;
    end
    if ~any(strcmp(metricName, innerMap.Properties.VariableNames))
        if any(strcmp('dV_DSM', innerMap.Properties.VariableNames))
            metricName = 'dV_DSM';
        else
            metricName = 'totalDv';
        end
    end

    x = rad2deg(innerMap.nu);
    y = innerMap.(metricName);
    if any(strcmp('isFeasible', innerMap.Properties.VariableNames))
        y(~innerMap.isFeasible) = NaN;
    end
    valid = isfinite(x) & isfinite(y);
    if ~any(valid)
        showEmptyAxes(ax, 'Inner search map has no finite values');
        hold(ax, 'off');
        return;
    end

    % Draw individual branch samples only. A line would falsely imply a
    % single-valued continuous curve through different branch families.
    h.points = scatter(ax, x(valid), y(valid), 24, y(valid), 'filled');
    colormap(ax, 'turbo');
    h.colorbar = colorbar(ax);
    ylabel(h.colorbar, metricName, 'Interpreter', 'none');
    xlabel(ax, '\nu_{DSM} [deg]');
    ylabel(ax, metricName, 'Interpreter', 'none');
    title(ax, 'Inner branch samples');
    xlim(ax, [0 360]);
    if strcmp(metricName, 'totalDv')
        ylim(ax, [0 100]);
    else
        ylim(ax, 'auto');
        xlim(ax, [0 360]);
    end

    marker = selectedMarker(selectedPoint, innerMap, metricName);
    if ~isempty(marker)
        h.selected(1) = plot(ax, marker(1), marker(2), 'ko', 'LineWidth', 1.8, 'MarkerSize', 12);
        h.selected(2) = plot(ax, marker(1), marker(2), 'rx', 'LineWidth', 2.2, 'MarkerSize', 12);
    end

    hold(ax, 'off');
end

function marker = selectedMarker(selectedPoint, innerMap, metricName)
    % When no explicit point is selected, mark the current best feasible
    % sample. Otherwise mark the point selected by click or auto-selection.
    marker = [];
    if isempty(selectedPoint)
        feasible = isfinite(innerMap.(metricName));
        if any(strcmp('isFeasible', innerMap.Properties.VariableNames))
            feasible = feasible & innerMap.isFeasible;
        end
        if ~any(feasible)
            return;
        end
        rows = innerMap(feasible, :);
        [~, idx] = min(rows.(metricName));
        selectedPoint = rows(idx, :);
    end

    if istable(selectedPoint)
        marker = [rad2deg(selectedPoint.nu(1)), selectedPoint.(metricName)(1)];
    elseif isstruct(selectedPoint) && isfield(selectedPoint, 'nu')
        marker = [rad2deg(selectedPoint.nu), selectedPoint.(metricName)];
    end
end

function showEmptyAxes(ax, message)
    text(ax, 0.5, 0.5, message, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'Color', [0.35 0.35 0.35]);
    xlabel(ax, '\nu_{DSM} [deg]');
    ylabel(ax, '\DeltaV [km/s]');
    title(ax, 'Inner branch samples');
    xlim(ax, [0 360]);
    ylim(ax, [0 100]);
end
